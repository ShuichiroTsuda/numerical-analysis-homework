FROM alpine

WORKDIR /home/fortran
RUN set -x && \
    apk update && \
    apk add --no-cache gfortran musl-dev

CMD ["/bin/sh"]