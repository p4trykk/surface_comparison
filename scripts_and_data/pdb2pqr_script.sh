#!/bin/bash

PDB_PATTERN="*_aligned.pdb"

echo "start"

for pdb in $PDB_PATTERN; do
    base=$(basename "$pdb" .pdb)

    cleaned="${base}_cleaned.pdb"
    top="${base}_cleaned.top"
    crd="${base}_cleaned.crd"
    pqr="${base}_cleaned.pqr"
    leap_in="${base}.in"

    echo "[1/4] pdb4amber"
    pdb4amber -i "$pdb" -o "$cleaned" 2>err.log

    if [ $? -ne 0 ]; then
        echo "  -> failed, trying --nohyd"
        pdb4amber -i "$pdb" -o "$cleaned" --nohyd 2>>err.log

        if [ $? -ne 0 ]; then
            echo "  -> -nohyd failed, trying --reduce"
            pdb4amber -i "$pdb" -o "$cleaned" --reduce
	    
            if [ $? -ne 0 ]; then
                echo "  -> pdb4amber FAILED, skipping $pdb"
                continue
            fi
        fi
    fi

    echo "[2/4] generating tleap input"

    cat > "$leap_in" <<EOF
source leaprc.protein.ff19SB
set default pbradii mbondi3
mol = loadpdb $cleaned
saveamberparm mol $top $crd
quit
EOF

    echo "[3/4] running tleap"
    tleap -f "$leap_in" > "${base}_tleap.log" 2>&1

    if [ ! -f "$top" ] || [ ! -f "$crd" ]; then
        echo "  -> tleap FAILED, skipping $pdb"
        continue
    fi

    echo "[4/4] generating PQR"
    ambpdb -p "$top" -c "$crd" -pqr > "$pqr"

    if [ $? -ne 0 ]; then
        echo "  -> ambpdb FAILED, skipping"
        continue
    fi

    echo "done $pqr"
done

echo "fin ;D"
