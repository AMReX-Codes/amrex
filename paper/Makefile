all:
	wget -N https://raw.githubusercontent.com/openjournals/whedon/master/resources/joss-logo.png
	wget -N https://raw.githubusercontent.com/openjournals/whedon/master/resources/latex.template
	pandoc \
    -V repository="https://github.com/amrex-codes/amrex" \
    -V archive_doi="http://dx.doi.org/10.5281/zenodo.2555438" \
    -V paper_url="http://joss.theoj.org/papers/10.xxxxx/joss.xxxxx" \
    -V formatted_doi="10.xxxxx/joss.xxxxx"  \
    -V review_issue_url="https://github.com/openjournals/joss-reviews/issues/588" \
    -V graphics="true" \
	-V journal_name="Journal of Open Source Software" \
    -V issue="34" \
    -V volume="4" \
    -V page="1265" \
    -V logo_path="joss-logo.png" \
    -V year="2018" \
    -V submitted="15 February 2019" \
    -V published="22 March 2019" \
    -V citation_author="Zhang et al." \
    -V paper_title="AMReX: a framework for block-structured adaptive mesh refinement (AMR)" \
    paper.md --template latex.template --filter pandoc-citeproc --bibliography paper.bib -o paper.pdf -V geometry:margin=1in

tex:
	wget -N https://raw.githubusercontent.com/openjournals/whedon/master/resources/joss-logo.png
	wget -N https://raw.githubusercontent.com/openjournals/whedon/master/resources/latex.template
	pandoc \
    -V repository="https://github.com/amrex-codes/amrex" \
    -V archive_doi="http://dx.doi.org/10.5281/zenodo.2555438" \
    -V paper_url="http://joss.theoj.org/papers/10.xxxxx/joss.xxxxx" \
    -V formatted_doi="10.xxxxx/joss.xxxxx" \
    -V review_issue_url="https://github.com/openjournals/joss-reviews/issues/588" \
    -V graphics="true" \
	-V journal_name="Journal of Open Source Software" \
    -V issue="34" \
    -V volume="4" \
    -V page="1265" \
    -V logo_path="joss-logo.png" \
    -V year="2018" \
    -V submitted="xx March 2019" \
    -V published="xx April 2019" \
    -V citation_author="Zhang et al." \
    -V paper_title="AMReX: a framework for block-structured adaptive mesh refinement (AMR)" \
    paper.md --template latex.template --filter pandoc-citeproc --bibliography paper.bib -o paper.tex -V geometry:margin=1in
