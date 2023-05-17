# Desing Project 2022-2023: Reference Bias

Contributors: Yunseol Park, Zhenyu Yang
Supervisors: Luca Renders, Lore Depyudt

`bwa_script.sh`: the complete alignment pipeline using BWA, `.out` file gives time, memory, and CPU usage.
`vg_script.sh`: the complete alignment pipeline using vg, `.out` file gives time, memory, and CPU usage.
`find_correctly_aligned.py`: caluclate rates of total aligned, correctly aligned, and incorrectly aligned reads out of the total number of reads.
`final_graph_all.py`: script for generating plots. The output of `find_correctly_aligned.py` is used as input.
