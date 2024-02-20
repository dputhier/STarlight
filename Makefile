MAKEFILE=Makefile
VERSION=0.0.1

.PHONY: help

#------------------------------------------------------------------
# Targets
#------------------------------------------------------------------


help:
	@echo ""
	@echo "- Available targets:"
	@echo "- Info: make check VERSION=0.0.1 "
	@perl -ne 'if(	/^([A-Za-z]+[A-Za-z_]*):/){print "\t",$$1,"\n"}' $(MAKEFILE)
	@echo ""
	@echo ""

clean:
	@rm -f src/*.o src/*.so; rm -f stcompr.Rcheck/dbfmcl/libs/dbfmcl.so; rm -rf ./dbfmcl.Rcheck; rm -rf ..Rcheck, rm -rf ./..Rcheck/
	@rm -rf /tmp/dbfmcl; rm -rf *dbf_out.txt; rm -rf *mcl_out.txt  rm -rf ./stcompr.Rcheck
	@rm -f tests/testthat/Rplot*; rm -rf tests/testthat/_snaps
	@rm -f *~

check: clean
	@rm -rf /tmp/stcompr; mkdir -p /tmp/stcompr; cp -r ./* /tmp/stcompr; cd /tmp/stcompr; \
	rm -f src/*.o src/*.so; rm -f stcompr.Rcheck/dbfmcl/libs/dbfmcl.so; \
	cd ..; R CMD build stcompr; R CMD check stcompr_$(VERSION).tar.gz

run_example:
	@echo "devtools::run_examples(pkg = '.')" | R --slave

checkfast: clean
	@rm -rf /tmp/stcompr; mkdir -p /tmp/stcompr; cp -r ./* /tmp/stcompr; cd /tmp/stcompr; \
	rm -f src/*.o src/*.so; rm -f stcompr.Rcheck/dbfmcl/libs/dbfmcl.so; \
	R CMD check --no-install .

doc:
	@echo ">>> Creating a package documentation"
	@echo "library(roxygen2); roxygen2::roxygenise()" | R --slave

install:
	@echo ">>> Installing..."
	@rm -f src/*.o src/*.so
	@R CMD INSTALL .

test:
	@echo ">>> Testing package"
	@rm -rf `ls tests/testthat/| grep -v \R$$`
	@echo "devtools::test()" | R --slave

test_by_file:
	@echo 'library(stcompr); for(i in dir("./tests/testthat/", pattern = ".R$$")){devtools::test_active_file(file.path("./tests/testthat/", i))}' | R --slave

coverage:
	@echo "Checking coverage"
	@echo "usethis::use_github_action('test-coverage'); cov <- covr::package_coverage(); print(as.data.frame(cov))" | R --slave

codecov:
	@echo "Uploading coverage (https://app.codecov.io/github/dputhier/stcompr)"
	@echo "library(covr); codecov(token ='8f08768a-0629-4ed0-91b9-bdd9f7019916')" | R --slave


#------------------------------------------------------------------
# Creating a release
#------------------------------------------------------------------

__check_defined_VER:
	@[ "$(VERSION)" ] || ( echo ">> VER is not set"; exit 1 )

release: __check_defined_VER
	@ echo "#-----------------------------------------------#"
	@ echo "# Starting the release $(VERSION)               #"
	@ echo "#-----------------------------------------------#"

release_bump: release
	@ echo "#-----------------------------------------------#"
	@ echo "# Bumping the program version                   #"
	@ echo "#-----------------------------------------------#"
	@ git checkout ./DESCRIPTION
	@ git checkout ./Makefile
	@ R CMD INSTALL .
	@ cat ./DESCRIPTION | perl -npe "s/Version: .*/Version: $(VERSION)/" > /tmp/stcompr.bump
	@ mv /tmp/stcompr.bump ./DESCRIPTION
	@ cat ./Makefile | perl -npe 's/^VERSION=.*/VERSION=$(VERSION)/' > /tmp/stcompr.bump
	@ mv /tmp/stcompr.bump ./Makefile
	@ echo "Version was bump to $(VERSION)"
	@ make install
	@ git commit -m 'Bumped version $(VERSION)'

readme:
	@ echo "- Rebuilting README.md from README.Rmd"
	@ echo "devtools::build_readme()" | R --slave

doc_html:
	@ echo "#-----------------------------------------------#"
	@ echo "# Building doc                                  #"
	@ echo "#-----------------------------------------------#"
	@ echo "Sys.setenv(RSTUDIO_PANDOC='/Applications/RStudio.app/Contents/Resources/app/quarto/bin/toolslibrary'); library(knitr); pkgdown::build_site()" | R --slave
	@ git add -u
	@ git commit -m "Updated html doc to $(VERSION)."


all: doc install check test



