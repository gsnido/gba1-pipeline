FROM mambaorg/micromamba:2.9.0

COPY --chown=$MAMBA_USER:$MAMBA_USER \
    environment.yml /tmp/environment.yml

RUN micromamba config set channel_priority strict \
    && micromamba install --yes --name base \
        --file /tmp/environment.yml \
    && micromamba clean --all --yes

COPY --chmod=0755 LONG-NEXT.sh /opt/gba1/bin/LONG-NEXT
COPY --chmod=0755 JOINT_GENOTYPING.sh /opt/gba1/bin/JOINT_GENOTYPING
COPY --chmod=0755 LONG-NEXT-SV.sh /opt/gba1/bin/LONG-NEXT-SV
COPY bedfiles/GBA1.bed /opt/gba1/GBA1.bed
COPY bedfiles/GBA1_GBA1LP_PSVs.hg38.tsv /opt/gba1/GBA1_GBA1LP_PSVs.hg38.tsv

ENV PATH="/opt/gba1/bin:${PATH}"


WORKDIR /work

CMD ["bash"]
