#!/bin/bash

# do some selection on date
DATE=`date --date="20 minutes ago" +"%Y-%m-%d %H:%M"`

docker run -ti --rm --env MOSAIC_DATE="$DATE" --env AWS_SECRET_ACCESS_KEY=$AWS_SECRET_ACCESS_KEY --env AWS_ACCESS_KEY_ID=$AWS_ACCESS_KEY_ID adokter/mosaic_aws Rscript /opt/mosaic_aws.R
