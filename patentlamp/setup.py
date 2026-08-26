from setuptools import setup, find_packages

setup(
    name="patentlamp",
    version="1.0.0",
    description="Autonomous Patent-Ready Molecular Engineering & Intellectual Property Engine for LAMP Primers",
    author="LaBiOmicS / NextLAMP Team",
    packages=find_packages(),
    entry_points={
        "console_scripts": [
            "patentlamp=patentlamp.cli:main",
        ],
    },
    python_requires=">=3.8",
)
