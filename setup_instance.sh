sudo yum update -y
sudo yum install docker -y
sudo service docker start
sudo usermod -a -G docker ec2-user
aws configure # provide keys
# log out and back in to get permissions right
