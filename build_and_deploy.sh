#!/bin/bash
set -e # Exit with nonzero exit code if anything fails

# First we test whether we can build libamrex or not
echo "Building libamrex... "
./configure
make -j4

# Then we build and deploy the sphinx / doxygen documentation
SOURCE_BRANCH="development"
TARGET_BRANCH="gh-pages"

# Pull requests and commits to other branches shouldn't try to deploy
if [ "$TRAVIS_PULL_REQUEST" != "false" -o "$TRAVIS_BRANCH" != "$SOURCE_BRANCH" ]; then
    echo "Skipping deploy."
    exit 0
fi

# Save some useful information
REPO=`git config remote.origin.url`
SSH_REPO=${REPO/https:\/\/github.com\//git@github.com:}
SHA=`git rev-parse --verify HEAD`

# Clone the existing gh-pages for this repo into out/
# Create a new empty branch if gh-pages doesn't exist yet (should only happen on first deploy)
git clone $REPO out
cd out
git checkout $TARGET_BRANCH || git checkout --orphan $TARGET_BRANCH
cd ..

# Clean out existing contents
rm -rf out/docs_html/**/* || exit 0
rm -rf out/tutorials_html/**/* || exit 0
rm -rf out/docs_xml/**/* || exit 0

# build the Doxygen documentation
cd Docs/Doxygen
doxygen doxygen.conf &> doxygen.out
cd ../..

mkdir -p out/docs_html
mkdir -p out/docs_xml
mkdir -p out/tutorials_html

# move it to the right place
mkdir -p out/docs_html/doxygen
mv Docs/Doxygen/html/* out/docs_html/doxygen/
mkdir -p out/docs_xml/doxygen
mv Docs/Doxygen/xml/* out/docs_xml/doxygen/

# run breathe and clean up
cd Docs/sphinx_documentation

#breathe-apidoc --o source ../../out/docs_xml/doxygen/ -g class,file
#python make_api.py

# now do sphinx
make SPHINX_BUILD="python -msphinx" latexpdf
mv build/latex/amrex.pdf source/ 
make SPHINX_BUILD="python -msphinx" html &> make_source_html.out

cd ../sphinx_tutorials
make SPHINX_BUILD="python -msphinx" latexpdf &> make_tutorials_latex.out
mv build/latex/amrex.pdf source/
make SPHINX_BUILD="python -msphinx" html &> make_tutorials_html.out
cd ../../

mv Docs/sphinx_documentation/build/html/* out/docs_html/
mv Docs/sphinx_tutorials/build/html/*     out/tutorials_html/
touch out/.nojekyll

# Now let's go have some fun with the cloned repo
cd out
git config user.name "Travis CI"
git config user.email "$COMMIT_AUTHOR_EMAIL"

if git diff-index --quiet HEAD; then
    exit 0
fi

# Commit the "changes", i.e. the new version.
# The delta will show diffs between new and old versions.
git add --all
git commit -m "Deploy to GitHub Pages: ${SHA}" || true

openssl aes-256-cbc -K $encrypted_6602cdd8f9c9_key -iv $encrypted_6602cdd8f9c9_iv -in ../id_rsa_travis.enc -out ../id_rsa_travis -d
chmod 600 ../id_rsa_travis
eval `ssh-agent -s`
ssh-add ../id_rsa_travis

git push $SSH_REPO $TARGET_BRANCH || true
ssh-agent -k
cd ..
