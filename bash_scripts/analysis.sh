

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_root="$(cd "$script_dir/.." && pwd)"

i=0
while read -r line; do
   line_arr=($line)
   range_i=${line_arr[0]} ind=$i envsubst < "$project_root/matplotlib_analysis/histogram.in" > "$project_root/tmp_cpptraj.in"
   $AMBERHOME/bin/cpptraj -i "$project_root/tmp_cpptraj.in" #>cpptraj.out
   
   ((i++))
# done < "$project_root/previous_iteration/iteration_3/cpptraj_base_pairing.txt"
done < "$project_root/Helix_separator/cpptraj_base_pairing.txt"
rm "$project_root/tmp_cpptraj.in"

xmgrace "$project_root/matplotlib_analysis/slide_probability_iteration_0.agr"