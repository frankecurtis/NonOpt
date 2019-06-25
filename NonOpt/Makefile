SUBDIRS = src problems exes tests

.PHONY: subdirs $(SUBDIRS)

subdirs: $(SUBDIRS)

$(SUBDIRS):
	$(MAKE) -C $@

clean:
	for d in $(SUBDIRS); \
	do \
	$(MAKE) --directory=$$d veryclean; \
	done

veryclean: clean
