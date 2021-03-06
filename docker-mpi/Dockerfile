FROM ubuntu:20.04
ARG BUILDER_UID
RUN DEBIAN_FRONTEND=noninteractive apt-get clean
RUN DEBIAN_FRONTEND=noninteractive apt-get update -y
RUN DEBIAN_FRONTEND=noninteractive apt-get install -y apt-utils vim bc build-essential pkgconf cmake m4 autoconf automake jq gdb python3 python3-dateutil apt-file curl less
RUN DEBIAN_FRONTEND=noninteractive apt-get install -y libmpich-dev
# ucspi-tcp: tcpserver used by tinysshd
RUN DEBIAN_FRONTEND=noninteractive apt-get install -y tinysshd ucspi-tcp  
# ssh client
RUN DEBIAN_FRONTEND=noninteractive apt-get install -y openssh-client
# ping, dig, nslookup, etc.
RUN DEBIAN_FRONTEND=noninteractive apt-get install -y iputils-ping dnsutils
RUN DEBIAN_FRONTEND=noninteractive apt-get install -y rsync

# Configure tinysshd
RUN tinysshd-makekey /etc/tinyssh/sshkeydir
RUN mkdir /root/.ssh/
COPY ssh/*  /root/.ssh/
RUN chmod 700 /root/.ssh/
RUN chmod 600 /root/.ssh/authorized_keys
RUN chmod 600 /root/.ssh/id_ecdsa

COPY bin-master/* /usr/local/bin/

# Embree
#
# Also you download Embree manually:
# wget https://github.com/embree/embree/releases/download/v3.12.1/embree-3.12.1.x86_64.linux.tar.gz
# and instead of ADD use this one:
#COPY embree-3.12.1.x86_64.linux.tar.gz /opt/
ADD https://github.com/embree/embree/releases/download/v3.12.1/embree-3.12.1.x86_64.linux.tar.gz /opt
RUN cd /opt && tar xzf embree-3.12.1.x86_64.linux.tar.gz
RUN echo "source /opt/embree-3.12.1.x86_64.linux/embree-vars.sh" >> /root/.profile

WORKDIR /root
CMD ["tcpserver", "-HRDl0", "0.0.0.0", "22", "/usr/sbin/tinysshd", "-v", "/etc/tinyssh/sshkeydir"]

