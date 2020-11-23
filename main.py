#!/usr/bin/env python3
#-*- coding:utf-8 -*-
u"""
Created at 2020.11.20
"""

import os

from glob import glob
from shutil import rmtree
from subprocess import check_call, CalledProcessError

import click
import pysam

from loguru import logger
from tqdm import tqdm


__dir__ = os.path.dirname(os.path.abspath(__file__))
__script__ = os.path.join(__dir__, "scripts")
__venv__ = os.path.join(__dir__, "venv")
PYTHON = os.path.join(__venv__, "bin/python")


class Config:

    def __init__(self, input_file, output, fasta, last_genome, ref_dict, name:str="", n_jobs:int=10):
        self.FASTQ = input_file if not input_file.endswith("bam") else ""
        self.BAM = input_file if input_file.endswith("bam") else ""
        if self.FASTQ and not name:
            self.name = os.path.basename(self.FASTQ).split(".")[0]
        elif self.BAM and not name:
            self.name = os.path.basename(self.BAM).split(".")[0]

        if name:
            self.name = name

        if not self.BAM:
            self.BAM = os.path.join(output, f"{name}.bam")
        elif not os.path.exists(self.BAM + ".bai"):
            logger.info("indexing %s" % self.BAM)
            pysam.index(self.BAM)

        self.CONSENSUS_FASTA = ""
        self.REFFASTA=fasta
        self.REFGENOME=last_genome
        self.REFDICT=ref_dict

        self.SV_CALLER = os.path.join(__venv__, "bin/NanoSV")
        self.SAMTOOLS="samtools"
        self.MINIMAP2=os.path.join(__venv__, "bin/minimap2")
        self.LAST_DIR=os.path.join(__venv__, "last-1145")
        self.WTDBG2_DIR=os.path.join(__venv__, "bin")
        self.PRIMER_DESIGN_DIR=os.path.join(__venv__, "primer3-2.5.0")

        self.n_jobs = n_jobs

        # output and temp file
        self.REGION_SELECTION_BED_OUTPUT = os.path.join(output, "regions.bed")
        self.REGION_SELECTION_BAM_OUTPUT = os.path.join(output, "regions.bam")

        self.VCF = os.path.join(output, f"{self.name}.vcf")
        self.VCF_FILTERED = ""
        self.CANDIDATE_DIR=os.path.join(output, "candidate_fusions")
        self.BAM_MERGE_OUT=os.path.join(output, "candidate_fusion_genes.bam")
        self.SV_CALLING_OUT=os.path.join(output, "candidate_fusion_genes.vcf")
        self.SV_CALLING_OUT_FILTERED = ""
        self.VCF_COMPLETE = os.path.join(output, "complete.vcf")

        self.FUSION_CHECK_VCF_OUTPUT = os.path.join(output, f"{self.name}_FusionGenes.vcf")
        self.FUSION_CHECK_INFO_OUTPUT=os.path.join(output, f"{self.name}_FusionGenesInfo.txt")
        self.FUSION_CHECK_PDF_OUTPUT=os.path.join(output, f"{self.name}_FusionGenes.pdf")

        self.PRIMER_DESIGN_BINDIR=os.path.join(self.PRIMER_DESIGN_DIR, "primers")
        self.PRIMER_DESIGN_GUIX_PROFILE=os.path.join(self.PRIMER_DESIGN_DIR, "emboss/.guix-profile")
        self.PRIMER_DESIGN_PRIMER3_CORE=os.path.join(self.PRIMER_DESIGN_DIR, "primer3/src/primer3_core")
        self.PRIMER_DESIGN_MISPRIMING=os.path.join(self.PRIMER_DESIGN_DIR, "repbase/current/empty.ref")
        self.PRIMER_DESIGN_PCR_TYPE='single'
        self.PRIMER_DESIGN_TILLING_PARAMS=''
        self.PRIMER_DESIGN_PSR='100-200'
        self.PRIMER_DESIGN_FLANK='200'

    def __str__(self):
        params = {
            "name": self.name,
            "FASTQ": self.FASTQ,
            "BAM": self.BAM,
            "REFFASTA": self.REFFASTA,
            "LAST_GENOME": self.REFGENOME,
            "REFDICT": self.REFDICT,
            "SV_CALLER": self.SV_CALLER,
            "SAMTOOLS": self.SAMTOOLS,
            "MINIMAP2": self.MINIMAP2,
            "LAST": self.LAST_DIR,
            "WTDBG2": self.WTDBG2_DIR,
            "PRIMER_DESIGN": self.PRIMER_DESIGN_DIR,
            "VCF": self.VCF
        }

        msg = [f"{k} - {v}" for k, v in params.items()]
        return "\n".join(msg)

    @property
    def temp(self):
        return [
            self.VCF_FILTERED,
            self.CANDIDATE_DIR,
            self.BAM_MERGE_OUT,
            self.BAM_MERGE_OUT + ".bai",
            self.SV_CALLING_OUT,
            self.SV_CALLING_OUT_FILTERED,
        ]

class Utils:

    @classmethod
    def __exists__(cls, path):
        if os.path.exists(path):
            if os.path.isfile(path):
                if os.path.getsize(path) > 0:
                    with open(path, "rb") as r:
                        try:
                            for rec in r:
                                pass
                            return True
                        except IOException as err:
                            pass
            else:
                if os.listdir(path):
                    return True
        return False

    @classmethod
    def call(cls, cmd: str, cwd=None, verbose:bool=True):
        if verbose:
            check_call(cmd, shell=True, cwd=cwd)
        else:
            with open(os.devnull) as w:
                check_call(cmd, shell=True, cwd=cwd, stdout=w, stderr=w)

    @classmethod
    def minimap_mapping(cls, config: Config, minimap_settings="-x map-ont -a --MD"):
        u"""
        call minimap2 mapping
        """
        logger.info("Mapping all reads using minimap2...")
        if config.FASTQ and os.path.exists(config.FASTQ):
            cls.call(f"{config.MINIMAP2} -t {config.n_jobs} {minimap_settings} {config.REF} {config.FASTQ} | {config.samtools} view -Sb -@ {config.n_jobs} | samtools sort -@ {config.n_jobs} > {config.BAM}")
            pysam.index(config.BAM)
        else:
            logger.info("Skip mapping")
        return config

    @classmethod
    def handle_selection(cls, config: Config, selection=None):

        if selection:
            logger.info("Seelcting regions to check for fusion genes...")

            try:
                cls.call(f"{PYTHON} {__script__}/RegionSelection.py -b {config.REGION_SELECTION_BED_OUTPUT} -r {selection}")
            except CalledProcessError as err:
                logger.error("!!! REGION SELECTION NOT CORRECTLY COMPLETED... exiting")
                exit(err)
        else:
            config.REGION_SELECTION_BAM_OUTPUT = config.BAM
            return config

        cls.call(f"{config.SAMTOOLS} view -H {config.BAM} ")

        with pysam.AlignmentFile(config.BAM) as r:
            with pysam.AlignmentFile(config.REGION_SELECTION_BAM_OUTPUT + ".unsorted", "wb+", template=r):
                with open(config.REGION_SELECTION_BED_OUTPUT) as r:
                    for line in r:
                        line = line.split()

                        for rec in r.fetch(config=line[0], start=int(line[1]), end=int(end[2])):
                            w.write(rec)

        pysam.sort("-o", config.REGION_SELECTION_BAM_OUTPUT, config.REGION_SELECTION_BAM_OUTPUT + ".unsorted")
        pysam.index(config.REGION_SELECTION_BAM_OUTPUT)

        if cls.__exists__(config.REGION_SELECTION_BAM_OUTPUT + ".unsorted"):
            os.remove(config.REGION_SELECTION_BAM_OUTPUT + ".unsorted")
        return config

    @classmethod
    def sv_calling(cls, config: Config, skip_first_calling: bool = False):
        logger.info("SV calling...")

        if not cls.__exists__(config.VCF) and not skip_first_calling:
            try:
                cls.call(f"{config.SV_CALLER} -s {config.SAMTOOLS} -c {os.path.join(__dir__, 'files/nanosv_last_config.ini')} -t {config.n_jobs} -o {config.VCF} {config.BAM}")
            except CalledProcessError as err:
                logger.error("failed to calling SV for %s" % config.VCF)
                exit(err)

        if cls.__exists__(config.VCF) and not config.VCF_FILTERED:
            config.VCF_FILTERED = cls.vcf_filter(config.VCF)

        if not cls.__exists__(config.SV_CALLING_OUT) and cls.__exists__(config.BAM_MERGE_OUT):
            logger.info("SV CALLING on merged BAM")
            try:
                cls.call(f"{config.SV_CALLER} -s {config.SAMTOOLS} -c {os.path.join(__dir__, 'files/nanosv_last_config.ini')} -t {config.n_jobs} -o {config.SV_CALLING_OUT} {config.BAM_MERGE_OUT}")
            except CalledProcessError as err:
                logger.error("failed to calling SV for %s" % config.SV_CALLING_OUT)
                exit(err)
        if cls.__exists__(config.SV_CALLING_OUT) and not config.SV_CALLING_OUT_FILTERED:
            config.SV_CALLING_OUT_FILTERED = cls.vcf_filter(config.SV_CALLING_OUT)

        return config

    @classmethod
    def vcf_filter(cls, input_vcf, filter:bool=True):
        if filter:
            logger.info("Removing insertions...")
            VCF_FILTERED = input_vcf.replace(".vcf", "_noINS.vcf")
            if not cls.__exists__(VCF_FILTERED):
                cls.call(f"grep \"^#\" {input_vcf} > {VCF_FILTERED}")
                cls.call(f"grep -v \"^#\" {input_vcf} | awk '$5!=\"<INS>\"' >> {VCF_FILTERED}")
        else:
            logger.info("Removing insertions and all variants without a PASS filter...")
            VCF_FILTERED = input_vcf.replace(".vcf", "_noINS_PASS.vcf")
            if not cls.__exists__(VCF_FILTERED):
                cls.call(f"grep \"^#\" {input_vcf} > {VCF_FILTERED}")
                cls.call(f"grep -v \"^#\" {input_vcf} >> {VCF_FILTERED}")
        return VCF_FILTERED

    @classmethod
    def extract_fusion(cls, config:Config, non_coding:bool):
        logger.info("Extracting reads that support candidate fusion genes...")

        if cls.__exists__(config.CANDIDATE_DIR):
            rmtree(config.CANDIDATE_DIR)

        os.makedirs(config.CANDIDATE_DIR, exist_ok=True)

        cmd = f"{PYTHON} {os.path.join(__script__, 'FusionReadExtraction.py')} -b {config.BAM} -v {config.VCF_FILTERED}  -o {config.CANDIDATE_DIR}"

        if non_coding:
            cmd = f"{cmd} -nc true"

        try:
            cls.call(cmd)
        except CalledProcessError as err:
            logger.error("!!! FUSION READ EXTRACTION NOT CORRECTLY COMPLETED... exiting")
            exit(err)

        return config

    @classmethod
    def consensus_calling(cls, config:Config, WTDBG2_SETTINGS='-x ont -g 3g', verbose: bool = False):

        if not config.CONSENSUS_FASTA:
            logger.info("Producing consensus if possible...")
            for fa in tqdm(glob(os.path.join(config.CANDIDATE_DIR, "*.fasta")), desc="consensus"):
                SVID = os.path.basename(fa).split("_")[0]
                prefix = fa.replace('.fasta', '_wtdbg2')

                try:
                    cls.call(f"{config.WTDBG2_DIR}/wtdbg2 {WTDBG2_SETTINGS} -i {fa} -t {config.n_jobs} -fo {prefix}", verbose=verbose)
                    cls.call(f"{config.WTDBG2_DIR}/wtpoa-cns -t {config.n_jobs} -i {prefix}.ctg.lay.gz -fo {prefix}.ctg.fa", verbose=verbose)
                except CalledProcessError as err:
                    logger.error(err)

                if cls.__exists__(f"{prefix}.ctg.fa"):
                    cls.call(f"sed -i \"s/>ctg/>${SVID}_ctg/g\" {prefix}.ctg.fa")
        else:
            logger.info("Skip consensus")
            consensus = os.path.join(config.WTDBG2_DIR, "consensus.fa")

            if not cls.__exists__(consensus):
                os.symlink(config.CONSENSUS_FASTA, consensus)
        return config

    @classmethod
    def mapping_fusion(cls, config:Config, use_last:bool, last_settings="", minimap_settings="-x map-ont -a --MD"):
        if use_last:
            logger.info("Using LAST for remapping")

            if not last_settings:
                last_settings = f'-Q 0 -p {config.LAST_DIR}/last_params'

            MAPPING_ARGS=f"-t {config.n_jobs} -r {config.REF} -rd {config.REFDICT} -l {config.LAST_DIR} -ls '' -s {config.SAMTOOLS}"

            for fa in tqdm(glob(os.path.join(config.WTDBG2_DIR, "*.fa")), desc="mapping fusion"):
                prefix = fa.replace(".fa", "")

                cls.call(f"{config.LAST_DIR}/src/lastal {last_settings} {config.REF} {fa} | {config.LAST_DIR}/src/last-split | {LAST_DIR}/scripts/maf-convert -f {config.REFDICT} sam /dev/stdin | {config.SAMTOOLS} view -b | {config.SAMTOOLS} sort > {prefix}.last.sorted.bam")
                pysam.index(f"{prefix}.last.sorted.bam")
        else:
            logger.info("Using minimap2 for remapping")

            for fa in glob(os.path.join(config.WTDBG2_DIR, "*.fa")):
                prefix = fa.replace(".fa", "")

                cls.call(f"{config.MINIMAP2} -t {config.n_jobs} {minimap_settings} {config.REF} {fa} | {config.samtools} view -Sb -@ {config.n_jobs} | samtools sort -@ {config.n_jobs} > {prefix}.sorted.bam")
                pysam.index(f"{prefix}.sorted.bam")

        fs = glob(os.path.join(config.WTDBG2_DIR, "*.sorted.bam"))

        if fs:
            template = pysam.AlignmentFile(fs[0])

            with pysam.AlignmentFile(config.BAM_MERGE_OUT + ".unsorted", "wb+", template=template) as w:
                for f in fs:
                    with pysam.AlignmentFile(f) as r:
                        for rec in r:
                            w.write(rec)

            template.close()

            pysam.sort("-o", config.BAM_MERGE_OUT, config.BAM_MERGE_OUT + ".unsorted")
            pysam.index(config.BAM_MERGE_OUT)

            if cls.__exists__(config.BAM_MERGE_OUT + ".unsorted"):
                os.remove(config.BAM_MERGE_OUT + ".unsorted")

        return config

    @classmethod
    def combine_sv(cls, config: Config):
        logger.info("Linking and combining SVs for complex fusion detection")

        if not cls.__exists__(config.VCF_COMPLETE):
            cls.call(f"{PYTHON} {__script__}/CombineSVs.py -v {config.SV_CALLING_OUT_FILTERED} -b {config.BAM_MERGE_OUT} -o {config.VCF_COMPLETE}")

        try:
            cls.call(f"grep -v \"^#\" {config.VCF_COMPLETE} >> {config.SV_CALLING_OUT_FILTERED}")
        except CalledProcessError as err:
            logger.warning(err)
        return config

    @classmethod
    def checking_fusion(cls, config: Config, non_coding:bool):
        try:
            cmd = f"{PYTHON} {__script__}/FusionCheck.py -ov {config.VCF} -o {config.FUSION_CHECK_VCF_OUTPUT} -fo {config.FUSION_CHECK_INFO_OUTPUT} -p {config.FUSION_CHECK_PDF_OUTPUT} -v {config.SV_CALLING_OUT_FILTERED} {'-nc' if non_coding else ''}"

            cls.call(cmd)
        except CalledProcessError as err:
            logger.error("!!! FUSION CHECK NOT CORRECTLY COMPLETED... exiting")
            exit(err)
        return config

    @classmethod
    def design_primers(cls, config: Config):
        primer_temp = os.path.join(config.PRIMER_DIR, "tmp")
        try:
            cmd = f"{PYTHON} {__script__}/PrimerFlankDesign.py -v {FUSION_CHECK_VCF_OUTPUT} -d {PRIMER_DIR} -f {PRIMER_DESIGN_FLANK}"

            cls.call(cmd)

            os.path.makedirs(primer_temp, exist_ok=True)

            for fa in glob(os.path.join(config.PRIMER_DIR, "*.fasta")):
                o = fa.replace(".fasta", ".primers")
                cls.call(f"bash {__dir__}/pipeline/primer_design.sh -f {fa} -o {o} -pdb {PRIMER_DESIGN_BINDIR} -pdpt {PRIMER_DESIGN_PCR_TYPE} -pdtp {PRIMER_DESIGN_TILLING_PARAMS} -psr {PRIMER_DESIGN_PSR} -pdgp {PRIMER_DESIGN_GUIX_PROFILE} -pdpc {PRIMER_DESIGN_PRIMER3_CORE} -pdm {PRIMER_DESIGN_MISPRIMING}", cwd=__dir__)

                out = os.path.join(primer_temp, "primer3.out")
                if cls.__exists__(out):
                    os.rename(out, fa.replace(".fasta", "_primerinfo.txt"))
        except CalledProcessError as err:
            logger.error("!!! PRIMER FLANK DESIGN NOT CORRECTLY COMPLETED... exiting")
            exit(err)
        finally:
            if cls.__exists__(primer_temp):
                rmtree(primer_temp)

    @classmethod
    def clean(cls, config: Config):
        for i in config.temp:
            if cls.__exists__(i):
                if os.path.isfile(i):
                    os.remove(i)
                elif os.path.isdir(i):
                    rmtree(i)


@click.command()
@click.option("-i", "--input-file", type=str, required=True, help="Path to fastq or bam file", show_default=True)
@click.option("-o", "--output", type=click.Path(), required=True, help="Path to output directory", show_default=True)
@click.option("-r", "--reference", type=click.Path(exists=True), required=True, help="Path to reference fasta file", show_default=True)
@click.option("-d", "--ref-dict", type=click.Path(), required=False, help="Path to reference dist file", show_default=True)
@click.option("-n", "--name", type=str, help="Output file name", show_default=True)
@click.option("-p", "--process", type=int, default=1, help="How many process to use", show_default=True)
@click.option("-s", "--selection", type=str, help="select specific region")
@click.option("--last-genome", type=click.Path(), required=False, help="Path to LAST reference genome, required by LAST", show_default=True)
@click.option("--use-last", is_flag=True, default=False, help="Wheter using LAST to align consensus sequence", show_default=True)
@click.option("--non-coding", is_flag=True, default=False, help="Also include non-coding fusions in the results (Not fully tested yet)", show_default=True)
@click.option("--not-filter", is_flag=True, default=False, help="Don't filter out all non-PASS SVs", show_default=True)
@click.option("--skip-first-sv", is_flag=True, default=False, help="If input bam is consensus then skip the processes before calling consensus", show_default=True)
def main(
    input_file: str, reference: str, ref_dict: str, last_genome: str,
    name: str, process: int, use_last: str, non_coding: bool, not_filter: bool,
    skip_first_sv: bool, selection: str, output: str
):
    if not ref_dict:
        ref_dict = reference + ".dict"
        Utils.call(f"gatk CreateSequenceDictionary -R {reference} -O {ref_dict}")

    config = Config(
        input_file=input_file,
        output=output,
        fasta=reference,
        ref_dict=ref_dict,
        last_genome=last_genome,
        name=name,
        n_jobs=process
    )
    logger.info(str(config))
    config = Utils.minimap_mapping(config)
    config = Utils.handle_selection(config, selection)

    if not skip_first_sv:
        config = Utils.sv_calling(config)
        config = Utils.extract_fusion(config, non_coding=non_coding)
        config = Utils.consensus_calling(config)
        config = Utils.mapping_fusion(config, use_last=use_last)
    else:
        config.BAM_MERGE_OUT = config.BAM
    config = Utils.sv_calling(config, skip_first_calling=skip_first_sv)
    config = Utils.combine_sv(config)
    config = Utils.checking_fusion(config, non_coding=non_coding)
    config = Utils.design_primers(config)


if __name__ == '__main__':
    main()
