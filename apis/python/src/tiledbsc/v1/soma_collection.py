from __future__ import annotations

from typing import Any, Dict, Iterator, Optional, Sequence

import tiledb

from .tiledb_array import TileDBArray
from .tiledb_object import TileDBObject
from .tiledb_platform_config import TileDBPlatformConfig
from .util import is_tiledb_creation_uri


class SOMACollection(TileDBObject):
    """
    Contains a key-value mapping where the keys are string names and the values are any SOMA-defined
    foundational or composed type, including SOMACollection, SOMADataFrame, SOMADenseNdArray,
    SOMASparseNdArray or SOMAExperiment.
    """

    # Cache to avoid repeated calls to the REST server for resolving group-member URIs
    # in the tiledb-cloud case. We invalidate this on add-member or remove-member.
    _cached_member_names_to_uris: Optional[Dict[str, str]]

    # TODO: comment re the performance impact of this cache.
    _cached_members: Dict[str, TileDBObject]

    def __init__(
        self,
        uri: str,
        *,
        name: Optional[str] = None,
        # Non-top-level objects can have a parent to propagate context, depth, etc.
        parent: Optional[SOMACollection] = None,
        # Top-level objects should specify these:
        tiledb_platform_config: Optional[TileDBPlatformConfig] = None,
        ctx: Optional[tiledb.Ctx] = None,
    ):
        """
        Also see the :class:`TileDBObject` constructor.
        """
        super().__init__(
            uri=uri,
            name=name,
            parent=parent,
            tiledb_platform_config=tiledb_platform_config,
            ctx=ctx,
        )
        self._cached_members = {}
        self._cached_member_names_to_uris = None

    def create(self) -> None:
        """
        Creates the data structure on disk/S3/cloud.
        """
        tiledb.group_create(uri=self._uri, ctx=self._ctx)
        self._common_create()  # object-type metadata etc

    # TODO
    # delete(uri)
    # Delete the SOMACollection specified with the URI.

    def __len__(self) -> int:
        """
        Returns the number of members in the collection.  Implements Python's `len(collection)`.
        """
        return len(self._get_member_names_to_uris())

    def __contains__(self, member_name: str) -> bool:
        """
        Tests for the existence of key in collection.
        Implements the `in` operator.
        """
        return member_name in self._get_member_names_to_uris()

    def get(self, member_name: str) -> TileDBObject:
        """
        Get the member object associated with the key
        """
        if member_name not in self._cached_members:
            # Do this here to avoid a cyclic module dependency:
            from .factory import _construct_member

            member_uri = self._get_child_uri(member_name)
            self._cached_members[member_name] = _construct_member(member_uri, self)

        return self._cached_members[member_name]

    def keys(self) -> Sequence[str]:
        """
        Gets the names of the members of the collection.
        """
        return self._get_member_names()

    def __getitem__(self, member_name: str) -> TileDBObject:
        """
        Get the member object associated with the key, when invoked as `collection["namegoeshere"]`.
        """
        return self.get(member_name)

    def __getattr__(self, member_name: str) -> TileDBObject:
        """
        Get the member object associated with the key, when invoked as `collection.namegoeshere`.
        """
        return self.get(member_name)

    def set(self, member: TileDBObject, *, relative: Optional[bool] = None) -> None:
        """
        Adds a member to the collection.
        """
        self._add_object(member, relative)
        self._cached_members[member.get_name()] = member

    def delete(self, member_name: str) -> None:
        """
        Removes a member from the collection, when invoked as `collection.delete("namegoeshere")`.
        """
        self._remove_object_by_name(member_name)

    def __delattr__(self, member_name: str) -> None:
        """
        Removes a member from the collection, when invoked as `del collection.namegoeshere`.
        """
        self.delete(member_name)

    def __delitem__(self, member_name: str) -> None:
        """
        Removes a member from the collection, when invoked as `del collection["namegoeshere"]`.
        """
        self.delete(member_name)

    def __iter__(self) -> Iterator[Any]:  # TODO: union type
        """
        Iterates over the collection.  Implements Python `for member in collection: ...` syntax.
        """
        for member_name, member_uri in self._get_member_names_to_uris().items():
            if member_name not in self._cached_members:
                # Do this here to avoid a cyclic module dependency:
                from .factory import _construct_member

                self._cached_members[member_name] = _construct_member(member_uri, self)
            yield self._cached_members[member_name]

    # ================================================================
    # PRIVATE METHODS FROM HERE ON DOWN
    # ================================================================

    def _tiledb_open(self, mode: str = "r") -> tiledb.Group:
        """
        This is just a convenience wrapper around tiledb group-open.
        It works asa `with self._tiledb_open() as G:` as well as `G = self._tiledb_open(); ...; G.close()`.
        """
        assert mode in ("r", "w")
        # This works in with-open-as contexts because tiledb.Group has __enter__ and __exit__ methods.
        return tiledb.Group(self._uri, mode=mode, ctx=self._ctx)

    def _create_unless_exists(self) -> None:
        """
        Auxiliary method for `_add_object`.
        """
        # Pre-checking if the group exists by calling tiledb.object_type is simple, however, for
        # tiledb-cloud URIs that occurs a penalty of two HTTP requests to the REST server, even
        # before a third, successful HTTP request for group-open.  Instead, we directly attempt the
        # group-create request, checking for an exception.
        try:
            self.create()
        except tiledb.cc.TileDBError as e:
            stre = str(e)
            # Local-disk/S3/tiledb-cloud exceptions all three say 'already exists'
            if "already exists" in stre:
                pass
            else:
                raise e

    def _get_child_uris(self, member_names: Sequence[str]) -> Dict[str, str]:
        """
        Computes the URIs for one or more children of the given object. For local disk, S3, and
        tiledb://.../s3://...  pre-creation URIs, this is simply the parent's URI, a slash, and the
        member name.  For post-creation TileDB-Cloud URIs, this is computed from the parent's
        information.  (This is because in TileDB Cloud, members have URIs like
        tiledb://namespace/df584345-28b7-45e5-abeb-043d409b1a97.)

        Since there's REST-server latency for getting name-to-URI mapping for group-member URIs, in
        the tiledb://... case, total latency is reduced when we ask for all group-element
        name-to-URI mappings in a single request to the REST server.
        """

        # TODO: TEMP
        if is_tiledb_creation_uri(self._uri):
            return self._get_child_creation_uris(member_names)

        # Pre-checking if the group exists by calling tiledb.object_type is simple, however, for
        # tiledb-cloud URIs that occurs a penalty of two HTTP requests to the REST server, even
        # before a third, successful HTTP request for group-open.  Instead, we directly attempt the
        # group-open request, checking for an exception.

        try:
            # If the group exists, get URIs for elements. In local disk, S3, etc this is just
            # concatenating a slash and the member name to the self URI; for TileDB Cloud,
            # the member URIs involve UUIDs which are known to the server.
            answer = {}

            mapping = self._get_member_names_to_uris()
            for member_name in member_names:
                if member_name in mapping:
                    answer[member_name] = mapping[member_name]
                else:
                    # Truly a slash, not os.path.join:
                    # * If the client is Linux/Un*x/Mac, it's the same of course
                    # * On Windows, os.path.sep is a backslash but backslashes are _not_ accepted for S3 or
                    #   tiledb-cloud URIs, whereas in Windows versions for years now forward slashes _are_
                    #   accepted for local-disk paths.
                    # This means forward slash is acceptable in all cases.
                    answer[member_name] = self._uri + "/" + member_name

            return answer

        except tiledb.cc.TileDBError as e:
            stre = str(e)
            # Local-disk/S3 does-not-exist exceptions say 'Group does not exist'; TileDB Cloud
            # does-not-exist exceptions are worded less clearly.
            if "Group does not exist" in stre or "HTTP code 401" in stre:
                return self._get_child_creation_uris(member_names)
            else:
                raise e

    def _get_child_creation_uris(self, member_names: Sequence[str]) -> Dict[str, str]:
        """
        Auxiliary method for `_get_child_uris`.
        """
        # Group being constructed. Here, appending "/" and name is appropriate in all cases:
        # even for tiledb://... URIs, pre-construction URIs are of the form
        # tiledb://namespace/s3://something/something/soma/membername.
        return {
            member_name: self._uri + "/" + member_name for member_name in member_names
        }

    def _get_child_uri(self, member_name: str) -> str:
        """
        Convenience wrapper around `_get_child_uris`.
        """
        return self._get_child_uris([member_name])[member_name]

    def _add_object(self, obj: TileDBObject, relative: Optional[bool] = None) -> None:
        """
        Adds a SOMA group/array to the current SOMA group -- e.g. base SOMA adding
        X, X adding a layer, obsm adding an element, etc.

        Semantics of `relative` from `self._tiledb_platform_config.member_uris_are_relative`:

        * If `False` then the group will have the absolute path of the member. For populating matrix
        elements within a SOMA in TileDB cloud, this is necessary. For populating SOMA elements within
        a SOMACollection on local disk, this can be useful if you want to be able to move the SOMACollection
        storage around and have it remember the (unmoved) locations of SOMA objects elsewhere, i.e.
        if the SOMACollectio is in one place while its members are in other places. If the SOMAs
        in the collection are contained within the SOMACollection directory, you probably want `relative=True`.

        * If `True` then the group will have the relative path of the member. For TileDB Cloud, this
        is never the right thing to do. For local-disk storage, this is essential if you want to move
        a SOMA to another directory and have it remember the locations of the members within it.

        * If `None`, then we select `relative=False` if the URI starts with `tiledb://`, else we
        select `relative=True`. This is the default.

        If the relative argument is supplied and is not None, it is used; secondly
        `self._tiledb_platform_config.member_uris_are_relative` is consulted; thirdly the URI prefix
        is consulted as described above.
        """

        if self._cached_member_names_to_uris is None:
            self._create_unless_exists()

        child_uri = obj.get_uri()
        child_name = obj.get_name()
        if relative is None:
            relative = self._tiledb_platform_config.member_uris_are_relative
        if relative is None:
            relative = not child_uri.startswith("tiledb://")
        if relative:
            child_uri = child_name
        self._cached_member_names_to_uris = None  # invalidate on add-member
        with self._tiledb_open("w") as G:
            G.add(uri=child_uri, relative=relative, name=child_name)
        # See _get_child_uri. Key point is that, on TileDB Cloud, URIs change from pre-creation to
        # post-creation. Example:
        # * Upload to pre-creation URI tiledb://namespace/s3://bucket/something/something/somaname
        # * Results will be at post-creation URI tiledb://namespace/somaname
        # * Note people can still use the pre-creation URI to read the data if they like.
        # * Member pre-creation URI tiledb://namespace/s3://bucket/something/something/somaname/obs
        # * Member post-creation URI tiledb://somaname/e4de581a-1353-4150-b1f4-6ed12548e497
        obj._uri = self._get_child_uri(child_name)

    def _remove_object(self, obj: TileDBObject) -> None:
        self._remove_object_by_name(obj.get_name())

    def _remove_object_by_name(self, member_name: str) -> None:
        self._cached_member_names_to_uris = None  # invalidate on remove-member
        if self._uri.startswith("tiledb://"):
            mapping = self._get_member_names_to_uris()
            if member_name not in mapping:
                raise Exception(f"name {member_name} not present in group {self._uri}")
            member_uri = mapping[member_name]
            with self._tiledb_open("w") as G:
                G.remove(member_uri)
        else:
            with self._tiledb_open("w") as G:
                G.remove(member_name)

    def _get_member_names(self) -> Sequence[str]:
        """
        Returns the names of the group elements. For a SOMACollection, these will SOMA names;
        for a SOMA, these will be matrix/group names; etc.
        """
        return list(self._get_member_names_to_uris().keys())

    def _get_member_uris(self) -> Sequence[str]:
        """
        Returns the URIs of the group elements. For a SOMACollection, these will SOMA URIs;
        for a SOMA, these will be matrix/group URIs; etc.
        """
        return list(self._get_member_names_to_uris().values())

    def _get_member_names_to_uris(self) -> Dict[str, str]:
        """
        Like `_get_member_names()` and `_get_member_uris`, but returns a dict mapping from
        member name to member URI.
        """
        if self._cached_member_names_to_uris is None:
            with self._tiledb_open("r") as G:
                self._cached_member_names_to_uris = {obj.name: obj.uri for obj in G}
        return self._cached_member_names_to_uris

    def _show_metadata(self, recursively: bool = True, indent: str = "") -> None:
        """
        Shows metadata for the group, recursively by default.
        """
        print(f"{indent}[{self._name}]")
        for key, value in self.metadata:  # XXX TEMP
            print(f"{indent}- {key}: {value}")
        if recursively:
            child_indent = indent + "  "
            with self._tiledb_open() as G:
                for obj in G:  # This returns a tiledb.object.Object
                    # It might appear simpler to have all this code within TileDBObject class,
                    # rather than (with a little duplication) in SOMACollection and TileDBArray.
                    # However, getting it to work with a recursive data structure and finding the
                    # required methods, it was simpler to split the logic this way.
                    object_type = tiledb.object_type(obj.uri, ctx=self._ctx)
                    if object_type == "group":
                        group = SOMACollection(uri=obj.uri, name=obj.name, parent=self)
                        group._show_metadata(recursively, indent=child_indent)
                    elif object_type == "array":
                        array = TileDBArray(uri=obj.uri, name=obj.name, parent=self)
                        array._show_metadata(recursively, indent=child_indent)
                    else:
                        raise Exception(
                            f"Unexpected object_type found: {object_type} at {obj.uri}"
                        )
