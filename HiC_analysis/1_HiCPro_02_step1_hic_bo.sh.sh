FASTQFILE=$SLURM_SUBMIT_DIR/inputfiles_hic_bo.txt; export FASTQFILE
make --file /share/lab_teng/trainee/JingjingWu/EBV/HiC_performance/HiC-Pro_install_path/HiC-Pro_3.1.0/scripts/Makefile CONFIG_FILE=/share/lab_teng/trainee/JingjingWu/EBV/hicpro_mapping_2/0_config-hicpro.txt CONFIG_SYS=/share/lab_teng/trainee/JingjingWu/EBV/HiC_performance/HiC-Pro_install_path/config-system.txt mapping  2>&1
~
