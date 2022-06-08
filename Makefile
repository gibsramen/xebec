TMPDIR := $(shell mktemp -d)
COOKIE_DIR := $(shell pwd)/xebec
TABLE_FILE := $(shell realpath xebec/tests/data/table.biom)
MD_FILE := $(shell realpath xebec/tests/data/metadata.tsv)
TREE_FILE := $(shell realpath xebec/tests/data/tree.tre)

all: test snaketest

test:
	@pytest xebec/tests

snaketest:
	@cd $(TMPDIR); \
	xebec -ft $(TABLE_FILE) -m $(MD_FILE) -t $(TREE_FILE) -o ./diversity-benchmark; \
	cd diversity-benchmark; \
	snakemake --cores 1; \
	snakemake --report report.html; \
	rm -r $(TMPDIR)
