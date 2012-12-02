#!/bin/bash
#
rm mesh_bandwidth_prb_output.txt
#
cp ../../data/polygonal_surface/sphere_q4_elements.txt sphere_q4_elements.txt
~/bincpp/$ARCH/mesh_bandwidth sphere_q4_elements.txt > mesh_bandwidth_prb_output.txt
rm sphere_q4_elements.txt
#
cp ../../data/polygonal_surface/sphere_t3_elements.txt sphere_t3_elements.txt
~/bincpp/$ARCH/mesh_bandwidth sphere_t3_elements.txt >> mesh_bandwidth_prb_output.txt
rm sphere_t3_elements.txt
#
cp ../../data/tet_mesh_order4/cube_elements.txt cube_order4_elements.txt
~/bincpp/$ARCH/mesh_bandwidth cube_order4_elements.txt >> mesh_bandwidth_prb_output.txt
rm cube_order4_elements.txt
#
cp ../../datasets/tet_mesh_order4/twenty_elements.txt twenty_order4_elements.txt
~/bincpp/$ARCH/mesh_bandwidth twenty_order4_elements.txt >> mesh_bandwidth_prb_output.txt
rm twenty_order4_elements.txt
#
cp ../../data/tet_mesh_order10/cube_order10_elements.txt cube_order10_elements.txt
~/bincpp/$ARCH/mesh_bandwidth cube_order10_elements.txt >> mesh_bandwidth_prb_output.txt
rm cube_order10_elements.txt
#
cp ../../data/tet_mesh_order10/oneoneeight_order10_elements.txt oneoneeight_order10_elements.txt
~/bincpp/$ARCH/mesh_bandwidth oneoneeight_order10_elements.txt >> mesh_bandwidth_prb_output.txt
rm oneoneeight_order10_elements.txt
#
cp ../../datasets/triangulation_order3/ell3_elements.txt ell3_elements.txt
~/bincpp/$ARCH/mesh_bandwidth ell3_elements.txt >> mesh_bandwidth_prb_output.txt
rm ell3_elements.txt
#
cp ../../datasets/triangulation_order3/hex_holes3_elements.txt hex_holes3_elements.txt
~/bincpp/$ARCH/mesh_bandwidth hex_holes3_elements.txt  >> mesh_bandwidth_prb_output.txt
rm hex_holes3_elements.txt
#
cp ../../datasets/triangulation_order3/hot_pipe3_elements.txt hot_pipe3_elements.txt
~/bincpp/$ARCH/mesh_bandwidth hot_pipe3_elements.txt >> mesh_bandwidth_prb_output.txt
rm hot_pipe3_elements.txt
#
cp ../../datasets/triangulation_order6/ell6_elements.txt ell6_elements.txt
~/bincpp/$ARCH/mesh_bandwidth ell6_elements.txt >> mesh_bandwidth_prb_output.txt
rm ell6_elements.txt
#
cp ../../datasets/triangulation_order6/hex_holes6_elements.txt hex_holes6_elements.txt
~/bincpp/$ARCH/mesh_bandwidth hex_holes6_elements.txt >> mesh_bandwidth_prb_output.txt
rm hex_holes6_elements.txt
#
echo "Program output written to mesh_bandwidth_prb_output.txt"
