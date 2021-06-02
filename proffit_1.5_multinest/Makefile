
# settings
CMAKE = cmake
#

BUILDDIR := $(CURDIR)/build


# default target
all:

# deep clean
purge:
	$(RM) -r $(BUILDDIR)
	find $(CURDIR) -name "*.pyc" -exec $(RM) -v \{} \;

# delegate any target to the build directory (except 'purge')
ifneq ($(MAKECMDGOALS),purge)
%: $(BUILDDIR)/Makefile FORCE
	$(MAKE) -C build $*
endif

# aliases
.PHONY: configure tests FORCE
ifneq ($(wildcard $(BUILDDIR)/Makefile),)
configure: rebuild_cache
else
configure: $(BUILDDIR)/Makefile
endif
	@ # do not delegate further

tests: all
	-$(MAKE) -C build test

# ensure that the target are always passed to the CMake Makefile
FORCE:
	@ # dummy target

# Special trick to allow a non-standard makefile name
#  If the makefile is not called 'Makefile', we get its update delegated to
#  cmake, unless we block the delegation.
$(lastword $(MAKEFILE_LIST)):
	@ # do not delegate further

# trigger CMake configuration
$(BUILDDIR)/Makefile:
	mkdir -p $(BUILDDIR)
	cd $(BUILDDIR) && $(CMAKE) $(CMAKEFLAGS) $(CURDIR)
