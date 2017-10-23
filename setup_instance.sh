sudo yum update -y
sudo yum install docker -y
sudo service docker start
sudo usermod -a -G docker ec2-user
aws configure # provide keys
# log out and back in to get permissions right
crontab -e
# for example:
# 0,10,20,30,40,50 0-13,22,23 * * * /home/ec2-user/mosaic.sh >/dev/null 2>&1
# note that a newline is required at the end of the statement
