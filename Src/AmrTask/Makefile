include ./arch.common

SUBDIRS = graph $(RTS_DIR) tutorials/UnitTests/

.PHONY: build $(SUBDIRS)

build: $(SUBDIRS)

$(SUBDIRS):
	$(MAKE) -C $@

.PHONY: all

all: build

clean:
	$(foreach dir, $(SUBDIRS), $(MAKE) -C $(dir) clean;)
