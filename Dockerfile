FROM uwgac/topmed-master:2.6.0

#Install METAL
RUN mkdir METAL && \
    cd METAL && \
    wget http://csg.sph.umich.edu/abecasis/Metal/download/Linux-metal.tar.gz && \
    gzip -d Linux-metal.tar.gz && \
    tar -xvf Linux-metal.tar && \
    mv generic-metal/metal /home/analyst && \
    cd .. && \
    rm -r METAL

# Install scripts
COPY *.sh /home/analyst/
RUN chmod a+x /home/analyst/*.sh
COPY *.R /home/analyst/
