from setuptools import setup, find_packages

setup(
    name="nextlamp",
    version="1.0.0",
    description="NextLAMP: A Modern, FAIR-Compliant, High-Performance Whole-Genome LAMP Primer Design Tool",
    author="Antigravity Team",
    packages=find_packages(include=["nextlamp", "nextlamp.*", "scripts", "scripts.*"]),
    package_data={
        "nextlamp": [
            "sample_data/*",
            "tests/glapd_comparison/*",
        ],
    },
    python_requires=">=3.10",
    install_requires=[
        "biopython>=1.80",
        "requests>=2.28.0",
    ],
    entry_points={
        "console_scripts": [
            "nextlamp=scripts.nextlamp.run_nextlamp:main",
            "nextlamp-prep=scripts.nextlamp.prep_data:main",
            "nextlamp-test=nextlamp.tests.test_subsample:main",
            "nextlamp-compare=nextlamp.tests.glapd_comparison.run_comparison:main",
            "nextlamp-compare-loops=nextlamp.tests.glapd_comparison.run_parallel_loop_evaluation:main",
        ],
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: POSIX :: Linux",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
)
