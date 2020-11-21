#!/usr/bin/env python
import os
from subprocess import check_call


__dir__ = os.path.dirname(os.path.abspath(__file__))


def build_minimap2():
    if not os.path.exists(os.path.join(__dir__, "venv/minimap2-2.17")):
        check_call(f"tar -C venv -xjf softwares/minimap2*.tar.bz2", shell=True, cwd=__dir__)
    
    if not os.path.exists(os.path.join(__dir__, "venv/minimap2-2.17/minimap2")):
        check_call("make", cwd=os.path.join(__dir__, "venv/minimap2-2.17"), shell=True)
    
    if not os.path.exists(os.path.join(__dir__, "venv/bin/minimap2")):
        os.symlink(os.path.join(__dir__, "venv/minimap2-2.17/minimap2"), os.path.join(__dir__, "venv/bin/minimap2"))


def build_wtdbg():
    check_call(f"tar -C venv -xzf softwares/wtdbg-2.5_x64_linux.tgz", shell=True, cwd=__dir__)
    check_call(f"mv venv/wtdbg-2.5_x64_linux/* venv/bin && rm -r venv/wtdbg-2.5_x64_linux/", shell=True, cwd=__dir__)

def build_last():
    if not os.path.exists(os.path.join(__dir__, "venv/last-1145")):
        check_call(f"unzip -d venv/  softwares/last-1145.zip", shell=True, cwd=__dir__)
        check_call("make", shell=True, cwd=os.path.join(__dir__, "venv/last-1145"))


def build_primer3():
    if not os.path.exists(os.path.join(__dir__, "venv/primer3-2.5.0")):
        check_call(f"tar -C venv -xzf softwares/primer3-v2.5.0.tar.gz", shell=True, cwd=__dir__)
        check_call("make", shell=True, cwd=os.path.join(__dir__, "venv/primer3-2.5.0/src"))


if __name__ == '__main__':
    if not os.path.exists(os.path.join(__dir__, "venv")):
        check_call(f"pip install -i https://mirrors.aliyun.com/pypi/simple virtualenv", shell=True)
        check_call(f"virtualenv venv -p python3", cwd=__dir__, shell=True)

    check_call(f"venv/bin/pip install -i https://mirrors.aliyun.com/pypi/simple certifi chardet matplotlib nltk numpy pysam PyVCF requests pybiomart NanoSV tqdm loguru click", shell=True)

    build_minimap2()
    build_wtdbg()
    build_last()
