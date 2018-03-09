sudo yum update -y
sudo yum install docker -y
sudo service docker start
sudo usermod -a -G docker ec2-user
aws configure # provide keys
# log out and back in to get permissions right
crontab -e
# for example:
# 5,15,25,35,45,55 0-16,22,23 * * * /home/ec2-user/mosaic.sh >/dev/null 2>&1
# note that a newline is required at the end of the crontab statement
