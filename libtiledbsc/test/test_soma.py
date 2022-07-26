import os
import pandas as pd
import pyarrow as pa
import libtiledbsc as sc

TEST_DIR = os.path.dirname(__file__)
SOMA_URI = f"{TEST_DIR}/../../test/soco/pbmc3k_processed"


def test_soma_list():
    soma = sc.SOMA(SOMA_URI)
    array_uris = soma.list_arrays()
    assert len(array_uris) == 19


def test_soma_query():
    sc.config_logging("debug", "soma.log")

    for i in range(22, 26):
        config = {"soma.init_buffer_bytes": f"{1 << i}"}
        sc.debug()
        sc.debug(f"SOMA query: {config}")
        sc.debug()
        soma = sc.SOMA(SOMA_URI, config)
        sq = soma.query()

        table = sq.next_results()
        while chunk := sq.next_results():
            table = pa.concat_tables([table, chunk])

        assert len(table.to_pandas()) == 4848644


if __name__ == "__main__":
    test_soma_query()
