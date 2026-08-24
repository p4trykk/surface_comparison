#!/bin/bash

RAD_VALUES=(0.5 1.0 1.5 2.0 2.5 3.0 3.5)

for file in *_4rev_aligned.pdb; do
    [ -e "$file" ] || { echo "cannot find files *_4rev_aligned.pdb"; exit 1; }
    base_name="${file%.pdb}"
    maps_dir="output_SURFMAP_${base_name}_bfactor/maps"
    echo "< $file >"
    for rad in "${RAD_VALUES[@]}"; do
        #echo "-rad $rad without norm"
        surfmap -pdb "$file" -tomap bfactor --png --nosmooth --keep --clean-map -rad "$rad" -s 2 --write-mapped-shell-pdb
        if [ -d "$maps_dir" ]; then
            mv "$maps_dir/${base_name}_map.png" "maps_vis_wnorm/${base_name}_map_rad_${rad}.png" 2>/dev/null
            mv "$maps_dir/${base_name}_map_clean.png" "without_norm/${base_name}_map_clean_rad_${rad}.png" 2>/dev/null
        else
            echo "paath $maps_dir not exist"
        fi

        #echo "-rad $rad  with norm"
        surfmap -pdb "$file" -tomap bfactor --png --nosmooth --keep --clean-map -rad "$rad" -s 2 --write-mapped-shell-pdb --bfactor-normalize

        if [ -d "$maps_dir" ]; then
            mv "$maps_dir/${base_name}_map.png" "maps_vis_norm/${base_name}_map_rad_${rad}.png" 2>/dev/null
            mv "$maps_dir/${base_name}_map_clean.png" "norm/${base_name}_map_clean_rad_${rad}.png" 2>/dev/null
        else
            echo "path $maps_dir not exist"
        fi
    done
done

echo "endd"
