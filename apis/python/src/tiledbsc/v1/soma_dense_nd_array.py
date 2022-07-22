from typing import Any, List, Optional, Tuple, Union

import numpy as np
import pyarrow as pa
import tiledb

from .tiledb_array import TileDBArray
from .tiledb_group import TileDBGroup
from .util import tiledb_type_from_arrow_type


# TODO: rethink parenting -- add a middle layer
class SOMADenseNdArray(TileDBArray):
    """
    Represents ``X`` and others.
    """

    _shape: Tuple[int]

    def __init__(
        self,
        uri: str,
        *,
        name: Optional[str] = None,
        parent: Optional[TileDBGroup] = None,
    ):
        """
        Also see the :class:`TileDBObject` constructor.
        """
        super().__init__(uri=uri, name=name, parent=parent)

    def create(
        self,
        type: pa.DataType,
        shape: Union[Tuple, List[int]],
    ) -> None:
        """
        Create a SOMADenseNdArray named with the URI.

        :param type: an Arrow type defining the type of each element in the array. If the type is
        unsupported, an error will be raised.

        :param shape: the length of each domain as a list, e.g., [100, 10]. All lengths must be in
        the uint64 range.
        """

        # checks on shape
        assert len(shape) > 0
        for e in shape:
            assert e > 0

        level = self._tiledb_platform_config.string_dim_zstd_level

        dims = []
        for e in shape:
            dim = tiledb.Dim(
                # Use tiledb default names like `__dim_0`
                domain=(0, e - 1),
                tile=min(e, 2048),  # TODO: PARAMETERIZE
                dtype=np.uint64,
                filters=[tiledb.ZstdFilter(level=level)],
            )
            dims.append(dim)
        dom = tiledb.Domain(dims, ctx=self._ctx)

        attrs = [
            tiledb.Attr(
                name="data",
                dtype=tiledb_type_from_arrow_type(type),
                filters=[tiledb.ZstdFilter()],
                ctx=self._ctx,
            )
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

        self._common_create()  # object-type metadata etc

    # TODO: static/class method?
    #    def delete(uri: str) -> None
    #        """
    #        Delete the SOMADenseNdArray specified with the URI.
    #        """

    # TODO: static/class method?
    #    def exists(uri: str) -> bool
    #        """
    #        Return true if object exists and is a SOMADenseNdArray.
    #        """

    #    def get_metadata():
    #        """
    #        Access the metadata as a mutable [`SOMAMetadataMapping`](#SOMAMetadataMapping)
    #        """

    # get_type() is inherited from TileDBObject

    def get_shape(self) -> Tuple[int]:
        """
        Return length of each dimension, always a list of length ``ndims``
        """
        # TODO: cache read
        # return self._shape
        raise Exception("TBD")

    def get_ndims(self) -> int:
        """
        Return number of index columns
        """
        # TODO: cache read
        raise Exception("TBD")

    #    def get_schema(self) -> Arrow.Schema:
    #        """
    #        Return data schema, in the form of an Arrow Schema
    #        """

    def get_is_sparse(self) -> bool:
        """
        Returns ``False``.
        """
        return False

    def read(
        slice: Any,  # TODO: scalar, slice/range, or list of any of the above
        # TODO: partitions: Optional[SOMAReadPartitions] = None,,
        # TODO: result_order: one of 'row-major' or 'column-major'
    ) -> Any:  # TODO: Iterator[DenseReadResult]
        """
        Read a user-specified subset of the object, and return as one or more Arrow.Tensor.

        :param slice: per-dimension slice, expressed as a scalar, a range, or a list of both.

        :param partitions: an optional [`SOMAReadPartitions`](#SOMAReadPartitions) hint to indicate
        how results should be organized.

        :param result_order: order of read results. Can be one of row-major or column-major.

        The `read` operation will return a language-specific iterator over one or more Arrow Tensor
        objects and information describing them, allowing the incremental processing of results larger
        than available memory. The actual iterator used is delegated to language-specific SOMA specs. The
        `DenseReadResult` should include:

        * The coordinates of the slice (e.g., origin, shape)
        * an Arrow.Tensor with the slice values
        """
        raise Exception("TBD")

    def write(
        self,
        coords: Union[Tuple, List[int]],
        values: pa.Tensor,
    ) -> None:
        """
        Write an Arrow.Tensor to the persistent object. As duplicate index values are not allowed, index
        values already present in the object are overwritten and new index values are added.

        :param coords: location at which to write the tensor

        :param values: an Arrow.Tensor containing values to be written. The type of elements in `values` must
        match the type of the SOMADenseNdArray.
        """

        with self._tiledb_open("w") as A:
            A[coords] = values.to_numpy()
