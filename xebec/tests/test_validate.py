import os

import pytest

from xebec.src import _validate as vd


def test_validate_table(data_paths, tmp_path):
    err_biom = os.path.join(tmp_path, "err.biom")
    with open(err_biom, "w") as f:
        f.write("kachow")

    with pytest.raises(ValueError) as exc_info:
        vd.validate_table(err_biom)
    exp_err_msg = "Table is empty!"
    assert str(exc_info.value) == exp_err_msg

    with pytest.raises(FileNotFoundError) as exc_info:
        vd.validate_table("NOT A FILE")
    exp_err_msg = "[Errno 2] No such file or directory: 'NOT A FILE'"
    assert str(exc_info.value) == exp_err_msg


def test_validate_metadata(data_paths, tmp_path):
    err_md = os.path.join(tmp_path, "err.tsv")
    with open(err_md, "w") as f:
        f.write("kerblam")

    with pytest.raises(ValueError) as exc_info:
        vd.validate_metadata(err_md)
    exp_err_msg = "Metadata is empty!"
    assert str(exc_info.value) == exp_err_msg

    with pytest.raises(FileNotFoundError) as exc_info:
        vd.validate_metadata("NOT A FILE")
    exp_err_msg = "[Errno 2] No such file or directory: 'NOT A FILE'"
    assert str(exc_info.value) == exp_err_msg


def test_validate_metadata(data_paths, tmp_path):
    with pytest.raises(FileNotFoundError) as exc_info:
        vd.validate_tree("NOT A FILE")
    exp_err_msg = "[Errno 2] No such file or directory: 'NOT A FILE'"
    assert str(exc_info.value) == exp_err_msg
