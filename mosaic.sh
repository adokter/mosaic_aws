#!/bin/bash

# do some selection on date
DATE=`date --date="20 minutes ago" +"%Y-%m-%d %H:%M"`

ID=`sed -n 2p ~/.aws/credentials | awk '{print $3}'`
KEY=`sed -n 3p ~/.aws/credentials | awk '{print $3}'`

echo $ID
docker run --rm --env MOSAIC_DATE="$DATE" --env AWS_SECRET_ACCESS_KEY="$KEY" --env AWS_ACCESS_KEY_ID="$ID" adokter/mosaic_aws Rscript /opt/mosaic_aws.R
