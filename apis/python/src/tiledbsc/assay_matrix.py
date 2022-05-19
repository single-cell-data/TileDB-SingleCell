import tiledb
from .soma_options import SOMAOptions
from .tiledb_array import TileDBArray
from .tiledb_group import TileDBGroup
from .tiledb_object import TileDBObject
import tiledbsc.util as util

import scipy
import numpy as np

from typing import Optional
import time


class AssayMatrix(TileDBArray):
    """
    Wraps a TileDB sparse array with two string dimensions.
    Used for X, obsp members, and varp members.
    """

    row_dim_name: str  # obs_id for X, obs_id_i for obsp; var_id_i for varp
    col_dim_name: str  # var_id for X, obs_id_j for obsp; var_id_j for varp
    attr_name: str

    # ----------------------------------------------------------------
    def __init__(
        self,
        uri: str,
        name: str,
        row_dim_name: str,
        col_dim_name: str,
        parent: Optional[TileDBGroup] = None,
    ):
        """
        See the TileDBObject constructor.
        """
        super().__init__(uri=uri, name=name, parent=parent)

        self.row_dim_name = row_dim_name
        self.col_dim_name = col_dim_name
        self.attr_name = "value"

    # ----------------------------------------------------------------
    # We don't have a .shape() method since X is sparse. One should
    # instead use the row-counts for the soma's obs and var.

    # ----------------------------------------------------------------
    def from_matrix(self, matrix, row_names, col_names) -> None:
        """
        Imports a matrix -- nominally scipy.sparse.csr_matrix or numpy.ndarray -- into a TileDB
        array which is used for X, obsp members, and varp members.
        """

        if self.verbose:
            s = util.get_start_stamp()
            print(f"{self.indent}START  WRITING {self.uri}")

        if self.exists():
            if self.verbose:
                print(f"{self.indent}Re-using existing array {self.uri}")
        else:
            self.create_empty_array(matrix_dtype=matrix.dtype)

        self.ingest_data(matrix, row_names, col_names)
        if self.verbose:
            print(util.format_elapsed(s, f"{self.indent}FINISH WRITING {self.uri}"))

    # ----------------------------------------------------------------
    def create_empty_array(self, matrix_dtype: np.dtype) -> None:
        """
        Create a TileDB 2D sparse array with string dimensions and a single attribute.
        """

        level = self.soma_options.string_dim_zstd_level
        dom = tiledb.Domain(
            tiledb.Dim(
                name=self.row_dim_name,
                domain=(None, None),
                dtype="ascii",
                filters=[tiledb.RleFilter()],
            ),
            tiledb.Dim(
                name=self.col_dim_name,
                domain=(None, None),
                dtype="ascii",
                filters=[tiledb.ZstdFilter(level=level)],
            ),
            ctx=self.ctx,
        )

        att = tiledb.Attr(
            self.attr_name,
            dtype=matrix_dtype,
            filters=[tiledb.ZstdFilter()],
            ctx=self.ctx,
        )

        sch = tiledb.ArraySchema(
            domain=dom,
            attrs=(att,),
            sparse=True,
            allows_duplicates=True,
            offsets_filters=[
                tiledb.DoubleDeltaFilter(),
                tiledb.BitWidthReductionFilter(),
                tiledb.ZstdFilter(),
            ],
            capacity=self.soma_options.X_capacity,
            cell_order=self.soma_options.X_cell_order,
            tile_order=self.soma_options.X_tile_order,
            ctx=self.ctx,
        )

        tiledb.Array.create(self.uri, sch, ctx=self.ctx)

    # ----------------------------------------------------------------
    def ingest_data(self, matrix, row_names, col_names) -> None:
        # TODO: add chunked support for CSC
        if (
            isinstance(matrix, scipy.sparse._csr.csr_matrix)
            and self.soma_options.write_X_chunked_if_csr
        ):
            self.ingest_data_rows_chunked(matrix, row_names, col_names)
        else:
            self.ingest_data_whole(matrix, row_names, col_names)

    # ----------------------------------------------------------------
    def ingest_data_whole(self, matrix, row_names, col_names) -> None:
        """
        Convert ndarray/(csr|csc)matrix to coo_matrix and ingest into TileDB.

        :param matrix: Matrix-like object coercible to a scipy coo_matrix.
        :param row_names: List of row names.
        :param col_names: List of column names.
        """

        assert len(row_names) == matrix.shape[0]
        assert len(col_names) == matrix.shape[1]

        mat_coo = scipy.sparse.coo_matrix(matrix)
        d0 = row_names[mat_coo.row]
        d1 = col_names[mat_coo.col]

        with tiledb.open(self.uri, mode="w", ctx=self.ctx) as A:
            A[d0, d1] = mat_coo.data

    # ----------------------------------------------------------------
    # Example: suppose this 4x3 is to be written in two chunks of two rows each
    # but written in sorted order.
    #
    # Original     Sorted     Permutation
    #  data       row names
    #
    #   X Y Z
    # C 0 1 2      A            1
    # A 4 0 5      B            2
    # B 7 0 0      C            0
    # D 0 8 9      D            3
    #
    # First chunk:
    # * Row indices 0,1 map to permutation indices 1,2
    # * i,i2 are 0,2
    # * chunk_coo is original matrix rows 1,2
    # * chunk_coo.row is [0,1]
    # * chunk_coo.row + i is [0,1]
    # * sorted_row_names: ['A', 'B']
    #
    # Second chunk:
    # * Row indices 2,3 map to permutation indices 0,3
    # * i,i2 are 2,4
    # * chunk_coo is original matrix rows 0,3
    # * chunk_coo.row is [0,1]
    # * chunk_coo.row + i is [2,3]
    # * sorted_row_names: ['C', 'D']
    #
    # See README-csr-ingest.md for important information of using this ingestor.
    # ----------------------------------------------------------------

    def ingest_data_rows_chunked(self, matrix, row_names, col_names) -> None:
        """
        Convert csr_matrix to coo_matrix chunkwise and ingest into TileDB.

        :param uri: TileDB URI of the array to be written.
        :param matrix: csr_matrix.
        :param row_names: List of row names.
        :param col_names: List of column names.
        """

        assert len(row_names) == matrix.shape[0]
        assert len(col_names) == matrix.shape[1]

        # Sort the row names so we can write chunks indexed by sorted string keys.  This will lead
        # to efficient TileDB fragments in the sparse array indexed by these string keys.
        #
        # Key note: only the _obs labels_ are being sorted, and along with them come permutation
        # indices for accessing the CSR matrix via cursor-indirection -- e.g. csr[28] is accessed as
        # with csr[permuation[28]] -- the CSR matrix itself isn't sorted in bulk.
        sorted_row_names, permutation = util.get_sort_and_permutation(list(row_names))
        # Using numpy we can index this with a list of indices, which a plain Python list doesn't support.
        sorted_row_names = np.asarray(sorted_row_names)

        s = util.get_start_stamp()
        if self.verbose:
            print(f"{self.indent}START  __ingest_coo_data_string_dims_rows_chunked")

        eta_tracker = util.ETATracker()
        with tiledb.open(self.uri, mode="w") as A:
            nrow = len(sorted_row_names)

            i = 0
            while i < nrow:
                t1 = time.time()
                # Find a number of CSR rows which will result in a desired nnz for the chunk.
                chunk_size = util.find_csr_chunk_size(
                    matrix, permutation, i, self.soma_options.goal_chunk_nnz
                )
                i2 = i + chunk_size

                # Convert the chunk to a COO matrix.
                chunk_coo = matrix[permutation[i:i2]].tocoo()

                # Write the chunk-COO to TileDB.
                d0 = sorted_row_names[chunk_coo.row + i]
                d1 = col_names[chunk_coo.col]

                if len(d0) == 0:
                    continue

                # Python ranges are (lo, hi) with lo inclusive and hi exclusive. But saying that
                # makes us look buggy if we say we're ingesting chunk 0:18 and then 18:32.
                # Instead, print doubly-inclusive lo..hi like 0..17 and 18..31.
                if self.verbose:
                    chunk_percent = 100 * (i2 - 1) / nrow
                    print(
                        "%sSTART  chunk rows %d..%d of %d (%.3f%%), obs_ids %s..%s, nnz=%d"
                        % (
                            self.indent,
                            i,
                            i2 - 1,
                            nrow,
                            chunk_percent,
                            d0[0],
                            d0[-1],
                            chunk_coo.nnz,
                        )
                    )

                # Write a TileDB fragment
                A[d0, d1] = chunk_coo.data

                if self.verbose:
                    t2 = time.time()
                    chunk_seconds = t2 - t1
                    eta = eta_tracker.ingest_and_predict(chunk_percent, chunk_seconds)

                    print(
                        "%sFINISH chunk in %.3f seconds, %7.3f%% done, ETA %s"
                        % (self.indent, chunk_seconds, chunk_percent, eta)
                    )

                i = i2

        if self.verbose:
            print(
                util.format_elapsed(
                    s, f"{self.indent}FINISH __ingest_coo_data_string_dims_rows_chunked"
                )
            )

    # ----------------------------------------------------------------
    def to_csr_matrix(self, row_labels, col_labels):
        """
        Reads the TileDB array storage for the storage and returns a sparse CSR matrix.  The
        row/columns labels should be obs,var labels if the AssayMatrix is X, or obs,obs labels if
        the AssayMatrix is obsp, or var,var labels if the AssayMatrix is varp.
        Note in all cases that TileDB will have sorted the row and column labels; they won't
        be in the same order as they were in any anndata object which was used to create the
        TileDB storage.
        """

        if self.verbose:
            s = util.get_start_stamp()
            print(f"{self.indent}START  read {self.uri}")

        # Since the TileDB array is sparse, with two string dimensions, we get back a dict:
        # * 'obs_id' key is a sequence of dim0 coordinates for X data.
        # * 'var_id' key is a sequence of dim1 coordinates for X data.
        # * 'values' key is a sequence of X data values.
        with tiledb.open(self.uri) as A:
            df = A[:]

        retval = util.X_and_ids_to_coo(
            df,
            self.row_dim_name,
            self.col_dim_name,
            self.attr_name,
            row_labels,
            col_labels,
        )

        if self.verbose:
            print(util.format_elapsed(s, f"{self.indent}FINISH read {self.uri}"))

        return retval
