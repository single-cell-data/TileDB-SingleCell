from typing import Any, Iterator, List, Optional, Sequence, TypeVar, Union

import numpy as np
import pyarrow as pa
import tiledb

from .tiledb_array import TileDBArray
from .tiledb_group import TileDBGroup
from .util import tiledb_type_from_arrow_type

ROWID = "soma_rowid"

Slice = TypeVar("Slice", bound=Sequence)


class SOMADataFrame(TileDBArray):
    """
    Represents ``obs``, ``var``, and others.

    A SOMADataFrame contains a "pseudo-column" called __rowid, of type uint64 and domain
    [0,num_rows).  The __rowid pseudo-column contains a unique value for each row in the
    SOMADataFrame, and is intended to act as a join key for other objects, such as a SOMANdArray.
    """

    # When we write to storage, we know if it's user-indexed or not. But when we are invoked
    # to read from storage, we don't know until we look at the storage.
    _is_user_indexed: Optional[bool]
    _index_column_names: Optional[List[str]]

    def __init__(
        self,
        uri: str,
        *,
        name: Optional[str] = None,
        parent: Optional[TileDBGroup] = None,
    ):
        """
        See also the :class:`TileDBOject` constructor.
        """
        super().__init__(uri=uri, name=name, parent=parent)
        self._is_user_indexed = None
        self._index_column_names = None

    def create(
        self,
        schema: pa.Schema,
        user_indexed: bool,
        index_column_names: Optional[List[str]] = None,
    ) -> None:
        """
        :param schema: Arrow Schema defining the per-column schema. This schema must define all
        columns, including columns to be named as index columns. The column name ``__rowid`` is
        reserved for the pseudo-column of the same name. If the schema includes types unsupported by
        the SOMA implementation, an error will be raised.

        :param user_indexed: If ``true``, creates a ``user-indexed`` dataframe, else a
        ``row-indexed`` dataframe.

        :param index_column_names: A list of column names to use as user-defined index columns
        (e.g., ``['cell_type', 'tissue_type']``). All named columns must exist in the schema, and at
        least one index column name is required. This parameter must be ``None`` if ``user_indexed``
        is ``False``.
        """
        if user_indexed:
            assert index_column_names is not None
            assert len(index_column_names) >= 1

            # assert all index_column_names are present in the schema
            schema_names_set = set(schema.names)
            for index_column_name in index_column_names:
                assert index_column_name in schema_names_set

            assert ROWID not in index_column_names
            assert ROWID not in schema_names_set

            self._create_empty_user_indexed(schema, index_column_names)
            self._is_user_indexed = True
            self._index_column_names = index_column_names
        else:
            assert index_column_names is None

            assert ROWID not in schema.names

            self._create_empty_row_indexed(schema)
            self._is_user_indexed = False
            self._index_column_names = []

        self._common_create()  # object-type metadata etc

    def _create_empty_user_indexed(
        self,
        schema: pa.Schema,
        index_column_names: Optional[List[str]],
    ) -> None:
        """
        Create a TileDB 1D sparse array with string dimension and multiple attributes.
        """

        level = self._tiledb_platform_config.string_dim_zstd_level

        dims = []
        for index_column_name in index_column_names:
            dtype = tiledb_type_from_arrow_type(schema.field(index_column_name).type)
            # We need domain=(None,None) for string dims
            lo = None
            hi = None
            if dtype != str:
                lo = np.iinfo(dtype).min
                hi = np.iinfo(dtype).max - 1
            dim = tiledb.Dim(
                name=index_column_name,
                domain=(lo, hi),
                tile=2048,  # TODO: PARAMETERIZE
                dtype=dtype,
                filters=[tiledb.ZstdFilter(level=level)],
            )
            dims.append(dim)

        dom = tiledb.Domain(dims, ctx=self._ctx)

        attrs = []

        attr = tiledb.Attr(
            name=ROWID,
            dtype=np.uint64,
            filters=[tiledb.ZstdFilter()],
            ctx=self._ctx,
        )
        attrs.append(attr)

        for attr_name in schema.names:
            if attr_name in index_column_names:
                continue
            attr = tiledb.Attr(
                name=attr_name,
                dtype=tiledb_type_from_arrow_type(schema.field(attr_name).type),
                filters=[tiledb.ZstdFilter()],
                ctx=self._ctx,
            )
            attrs.append(attr)

        sch = tiledb.ArraySchema(
            domain=dom,
            attrs=attrs,
            sparse=True,
            allows_duplicates=self._tiledb_platform_config.allows_duplicates,
            offsets_filters=[
                tiledb.DoubleDeltaFilter(),
                tiledb.BitWidthReductionFilter(),
                tiledb.ZstdFilter(),
            ],
            capacity=100000,
            cell_order="row-major",
            # As of TileDB core 2.8.2, we cannot consolidate string-indexed sparse arrays with
            # col-major tile order: so we write `X` with row-major tile order.
            tile_order="row-major",
            ctx=self._ctx,
        )

        tiledb.Array.create(self._uri, sch, ctx=self._ctx)

    def _create_empty_row_indexed(
        self,
        schema: pa.Schema,
    ) -> None:
        """
        Create a TileDB 1D dense array with uint64 __rowid dimension and multiple attributes.
        """

        level = self._tiledb_platform_config.string_dim_zstd_level

        dom = tiledb.Domain(
            tiledb.Dim(
                name=ROWID,
                domain=(0, np.iinfo(np.int64).max),
                tile=2048,  # TODO: PARAMETERIZE
                dtype=np.uint64,
                filters=[tiledb.ZstdFilter(level=level)],
            ),
            ctx=self._ctx,
        )

        attrs = [
            tiledb.Attr(
                name=attr_name,
                dtype=tiledb_type_from_arrow_type(schema.field(attr_name).type),
                filters=[tiledb.ZstdFilter()],
                ctx=self._ctx,
            )
            for attr_name in schema.names
        ]

        sch = tiledb.ArraySchema(
            domain=dom,
            attrs=attrs,
            sparse=False,
            allows_duplicates=self._tiledb_platform_config.allows_duplicates,
            offsets_filters=[
                tiledb.DoubleDeltaFilter(),
                tiledb.BitWidthReductionFilter(),
                tiledb.ZstdFilter(),
            ],
            capacity=100000,
            cell_order="row-major",
            tile_order="row-major",
            ctx=self._ctx,
        )

        tiledb.Array.create(self._uri, sch, ctx=self._ctx)

    # ----------------------------------------------------------------
    # TODO: static/class method?
    #    def delete(uri: str) -> None
    #        """
    #        Delete the SOMADataFrame specified with the URI.
    #        """

    # TODO: static/class method?
    #    def exists(uri: str) -> bool
    #        """
    #        Return true if object exists and is a SOMADataFrame.
    #        """

    #    def get_metadata():
    #        """
    #        Access the metadata as a mutable [`SOMAMetadataMapping`](#SOMAMetadataMapping)
    #        """

    # get_type() is inherited from TileDBObject

    #    def get_shape() -> Tuple[int]:
    #        """
    #        Return length of each dimension, always a list of length ``ndims``
    #        """

    def get_ndims(self) -> int:
        """
        Return number of index columns
        """
        return len(self.get_index_column_names())

    # TODO
    #    def get_schema(self) -> Arrow.Schema:
    #        """
    #        Return data schema, in the form of an Arrow Schema
    #        """

    def get_is_row_indexed(self) -> bool:
        """
        Return true if row-indexed, false if user-indexed.
        """
        return not self.get_is_user_indexed()

    def get_is_user_indexed(self) -> bool:
        """
        Return true if user-indexed, false if row-indexed.
        """
        # If we've cached the answer, skip the storage read. Especially if the storage is on the
        # cloud, where we'll avoid an HTTP request.
        if self._is_user_indexed is None:
            # TODO have written, and then read back, object-type metadata?
            # For now: we use sparse for user-indexed and dense for row-indexed,
            # so this check is a litmus test.
            #
            # If the array hasn't had its schema created then this will throw (as it should)
            # -- we have no way to answer the question with a true or a false.
            with self._tiledb_open() as A:
                self._is_user_indexed = A.schema.sparse
        return self._is_user_indexed

    def get_index_column_names(self) -> List[str]:
        """
        Return index (dimension) column names if user-indexed, or an empty list if row-indexed.
        """
        # If we've cached the answer, skip the storage read. Especially if the storage is on the
        # cloud, where we'll avoid an HTTP request.
        if self._index_column_names is None:
            if self.get_is_user_indexed():
                names = []
                with self._tiledb_open() as A:
                    dom = A.domain
                    for i in range(dom.ndim):
                        names.append(dom.dim(i).name)
                self._index_column_names = names
            else:
                self._index_column_names = []
        return self._index_column_names

    def read(
        self,
        *,
        # TODO: find out how to spell this in a way the type-checker will accept :(
        # ids: Optional[Union[Sequence[int], str, Slice]] = "all",
        ids: Optional[Any] = "all",
        column_names: Optional[Union[Sequence[str], str]] = "all",
        # TODO: partitions,
        # TODO: result_order,
        # TODO: value_filter
    ) -> Iterator[pa.RecordBatch]:
        """
        Read a user-defined subset of data, addressed by the dataframe indexing columns, optionally
        filtered, and return results as one or more Arrow.RecordBatch.

        :param ids: for each index dimension, which rows to read. Defaults to 'all'.

        :param column_names: the named columns to read and return. Defaults to 'all'.

        :param partitions: an optional ``SOMAReadPartitions`` hint to indicate how results should be
        organized.

        :param result_order: order of read results. If the dataframe is `user-indexed`, this can be one of
        'row-major', 'col-major', or 'unordered'. If the dataframe is 'row-indexed', this can be one of
        'rowid-ordered' or 'unordered'.

        :param value_filter: an optional [value filter] to apply to the results. Defaults to no
        filter.

        **Indexing**: the `ids` parameter will support, per dimension:

        * For `row-indexed` dataframes, a row offset (uint), a row-offset range (slice), or a list of both.
        * For `user-indexed` dataframes, a list of values of the type of the indexed column.
        """

        use_all_ids = False
        use_all_column_names = False
        if isinstance(ids, str):
            assert ids == "all"  # Enforce runtime type check
            use_all_ids = True
        if isinstance(column_names, str):
            assert column_names == "all"  # Enforce runtime type check
            use_all_column_names = True

        with self._tiledb_open("r") as A:
            q = A.query(return_arrow=True, return_incomplete=True)

            if use_all_ids:
                if use_all_column_names:
                    iterator = q.df[:]
                else:
                    iterator = q.df[:][column_names]
            else:
                if use_all_column_names:
                    iterator = q.df[ids]
                else:
                    iterator = q.df[ids][column_names]

            for df in iterator:
                batches = df.to_batches()
                for batch in batches:
                    yield batch

    def write(self, values: pa.RecordBatch) -> None:
        """
        Write an Arrow.RecordBatch to the persistent object. As duplicate index values are not allowed,
        index values already present in the object are overwritten and new index values are added.

        :param values: An Arrow.RecordBatch containing all columns, including the index columns. The
        schema for the values must match the schema for the :class:`SOMADataFrame`.

        If the dataframe is row-indexed, the ``values`` Arrow RecordBatch must contain a ``__rowid``
        (uint64) column, indicating which rows are being written.
        """
        if self.get_is_user_indexed():
            dim_cols_list = []
            attr_cols_map = {}
            dim_names_set = self.get_index_column_names()
            n = None

            for name in values.schema.names:
                n = len(values.column(name))
                if name in dim_names_set:
                    dim_cols_list.append(values.column(name).to_pandas())
                else:
                    attr_cols_map[name] = values.column(name).to_pandas()
            assert n is not None

            # XXX TEMP STUB
            attr_cols_map[ROWID] = np.asarray(range(n))

            dim_cols_list = [list(dim_col) for dim_col in dim_cols_list]
            with self._tiledb_open("w") as A:
                # TODO: find the right syntax for vardims ... it's not the `*` operator ...
                # A[*dim_cols_list] = attr_cols_map
                if len(dim_cols_list) == 1:
                    A[dim_cols_list[0]] = attr_cols_map
                elif len(dim_cols_list) == 2:
                    A[dim_cols_list[0], dim_cols_list[1]] = attr_cols_map
                else:
                    raise Exception("ndims >= 2 not currently supported")

        else:  # row-indexed

            assert ROWID in values.schema.names

            # TODO: contiguity check, and/or split into multiple contiguous writes
            # For now: just assert that these _already are_ contiguous and start with 0.
            rowids = [e.as_py() for e in values.column(ROWID)]
            assert len(rowids) > 0

            rowids = sorted(rowids)
            assert rowids[0] == 0

            lo = rowids[0]
            hi = rowids[-1]

            attr_cols_map = {}
            for name in values.schema.names:
                if name != ROWID:
                    attr_cols_map[name] = np.asarray(
                        values.column(name).to_pandas(
                            types_mapper=tiledb_type_from_arrow_type,
                        )
                    )

            with self._tiledb_open("w") as A:
                A[lo : (hi + 1)] = attr_cols_map
