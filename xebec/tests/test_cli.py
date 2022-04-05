from click.testing import CliRunner

from xebec.cli.cli import xebec


def test_xebec(data_paths):
    runner = CliRunner()
    with runner.isolated_filesystem():
        result = runner.invoke(xebec, [
            "-ft", data_paths.table_file,
            "-m", data_paths.metadata_file,
            "-t", data_paths.tree_file,
            "-o", "./test_output"
        ])
    assert result.exit_code == 0


def test_xebec_bad_validation(data_paths):
    runner = CliRunner()
    with runner.isolated_filesystem():
        result = runner.invoke(xebec, [
            "-ft", data_paths.table_file,
            "-m", "BAD_PATH",
            "-t", data_paths.tree_file,
            "-o", "./test_output",
        ])
    err_msg = str(result.exc_info[1])
    assert result.exit_code == 1
    assert isinstance(result.exception, FileNotFoundError)
    assert "No such file or directory" in err_msg


def test_xebec_no_validation(data_paths):
    runner = CliRunner()
    with runner.isolated_filesystem():
        result = runner.invoke(xebec, [
            "-ft", data_paths.table_file,
            "-m", "BAD_PATH",
            "-t", data_paths.tree_file,
            "--no-validate-input",
            "-o", "./test_output",
        ])
    assert result.exit_code == 0
