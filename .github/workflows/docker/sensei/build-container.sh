image_name=senseiinsitu/ci:fedora35-amrex-$(date +%Y%m%d)
docker build --tag $image_name .
docker push $image_name
