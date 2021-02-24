echo "start block TestFastaInput"

echo "  -testP3FastaInput"
java -jar target/primer3-1.0-SNAPSHOT.jar src/test/resources/fasta_input_test/primer_task_input -o test_results/testP3FastaInput_out -fasta src/test/resources/fasta_input_test/primer_fasta_input.fasta

diff test_results/testP3FastaInput_out src/test/resources/fasta_input_test/primer_task_output > test_results/diff_results/testP3FastaInput_diff

echo "  -testP3NInput"
java -jar target/primer3-1.0-SNAPSHOT.jar src/test/resources/fasta_input_test/p3_3_prime_n_input -o test_results/testP3NInput_out -fasta src/test/resources/fasta_input_test/p3_3_prime_n_input.fasta

diff test_results/testP3NInput_out src/test/resources/fasta_input_test/p3_3_prime_n_output > test_results/diff_results/testP3NInput_diff

echo "  -testP30Input"
java -jar target/primer3-1.0-SNAPSHOT.jar src/test/resources/fasta_input_test/p3_3_prime_0_input -o test_results/testP30Input_out -fasta src/test/resources/fasta_input_test/p3_3_prime_0_input.fasta

diff test_results/testP30Input_out src/test/resources/fasta_input_test/p3_3_prime_0_output > test_results/diff_results/testP30Input_diff

echo "end block TestFastaInput"
echo "start block TestGroupA"

echo "  -testP3ThermodAlignInput"
java -jar target/primer3-1.0-SNAPSHOT.jar src/test/resources/groupA_test/primer_thermod_align_input -o test_results/testP3ThermodAlignInput_out

diff test_results/testP3ThermodAlignInput_out src/test/resources/groupA_test/primer_thermod_align_output > test_results/diff_results/testP3ThermodAlignInput_diff

echo "  -testP3ThreePrimeDistanceInput"
java -jar target/primer3-1.0-SNAPSHOT.jar src/test/resources/groupA_test/primer_three_prime_distance_input -o test_results/testP3ThreePrimeDistanceInput_out

diff test_results/testP3ThreePrimeDistanceInput_out src/test/resources/groupA_test/primer_three_prime_distance_output > test_results/diff_results/testP3ThreePrimeDistanceInput_diff

echo "end block TestGroupA"
echo "start block TestGroupB"

echo "  -testP3Mispriming"
java -jar target/primer3-1.0-SNAPSHOT.jar src/test/resources/groupB_test/primer_mispriming_input -o test_results/testP3Mispriming_out

diff test_results/testP3Mispriming_out src/test/resources/groupB_test/primer_mispriming_output > test_results/diff_results/testP3Mispriming_diff

echo "  -testP3MisprimingTH"
java -jar target/primer3-1.0-SNAPSHOT.jar src/test/resources/groupB_test/primer_mispriming_th_input -o test_results/testP3MisprimingTH_out

diff test_results/testP3MisprimingTH_out src/test/resources/groupB_test/primer_mispriming_th_output > test_results/diff_results/testP3MisprimingTH_diff

echo "  -testP3MisprimingBoundary1"
java -jar target/primer3-1.0-SNAPSHOT.jar src/test/resources/groupB_test/primer_mispriming_boundary1_input -o test_results/testP3MisprimingBoundary1_out

diff test_results/testP3MisprimingBoundary1_out src/test/resources/groupB_test/primer_mispriming_boundary1_output > test_results/diff_results/testP3MisprimingBoundary1_diff

echo "  -testP3MisprimingBoundary2"
java -jar target/primer3-1.0-SNAPSHOT.jar src/test/resources/groupB_test/primer_mispriming_boundary2_input -o test_results/testP3MisprimingBoundary2_out

diff test_results/testP3MisprimingBoundary2_out src/test/resources/groupB_test/primer_mispriming_boundary2_output > test_results/diff_results/testP3MisprimingBoundary2_diff

echo "  -testP3MisprimingLongLib"
java -jar target/primer3-1.0-SNAPSHOT.jar src/test/resources/groupB_test/primer_mispriming_long_lib_input -o test_results/testP3MisprimingLongLib_out

diff test_results/testP3MisprimingLongLib_out src/test/resources/groupB_test/primer_mispriming_long_lib_output > test_results/diff_results/testP3MisprimingLongLib_diff

echo "  -testP3ObjFnInput"
java -jar target/primer3-1.0-SNAPSHOT.jar src/test/resources/groupB_test/primer_obj_fn_input -o test_results/testP3ObjFnInput_out

diff test_results/testP3ObjFnInput_out src/test/resources/groupB_test/primer_obj_fn_output > test_results/diff_results/testP3ObjFnInput_diff

echo "  -testP3LibAmbCodesInput"
echo "Ignored: Not Implemented"
# java -jar target/primer3-1.0-SNAPSHOT.jar src/test/resources/groupB_test/primer_lib_amb_codes_input -o test_results/testP3LibAmbCodesInput_out

#diff test_results/testP3LibAmbCodesInput_out src/test/resources/groupB_test/primer_lib_amb_codes_output > test_results/diff_results/testP3LibAmbCodesInput_diff

echo "end block TestGroupB"
echo "start block TestGroupC"

echo "  -testP3MaskerInput"
java -jar target/primer3-1.0-SNAPSHOT.jar src/test/resources/groupC_test/primer_masker_input -o test_results/testP3MaskerInput_out

diff test_results/testP3MaskerInput_out src/test/resources/groupC_test/primer_masker_output > test_results/diff_results/testP3MaskerInput_diff

echo "end block TestGroupC"
echo "start block TestMultiplexInput"

echo "  -testP3SpecificInput"
java -jar target/primer3-1.0-SNAPSHOT.jar src/test/resources/multiplex_test/primer_task_specific_input -o test_results/testP3SpecificInput_out -fasta src/test/resources/multiplex_test/primer_task_specific.fasta

diff test_results/testP3SpecificInput_out src/test/resources/multiplex_test/primer_task_specific_output > test_results/diff_results/testP3SpecificInput_diff

echo "  -testP3SpecificInput5T"
java -jar target/primer3-1.0-SNAPSHOT.jar src/test/resources/multiplex_test/primer_task_specific_input_5T -o test_results/testP3SpecificInput5T_out -fasta src/test/resources/multiplex_test/targets_new.fasta

diff test_results/testP3SpecificInput5T_out src/test/resources/multiplex_test/primer_task_specific_output > test_results/diff_results/testP3SpecificInput5T_diff

echo "  -testP3SpecificInput5Tm"
java -jar target/primer3-1.0-SNAPSHOT.jar src/test/resources/multiplex_test/primer_task_specific_input_5T_m -o test_results/testP3SpecificInput5Tm_out -fasta src/test/resources/multiplex_test/targets_new.fasta

diff test_results/testP3SpecificInput5Tm_out src/test/resources/multiplex_test/primer_task_specific_output > test_results/diff_results/testP3SpecificInput5Tm_diff

echo "  -testP3SpecificInput7T"
java -jar target/primer3-1.0-SNAPSHOT.jar src/test/resources/multiplex_test/primer_task_specific_input_7T -o test_results/testP3SpecificInput7T_out -fasta src/test/resources/multiplex_test/targets_new_fake.fasta

diff test_results/testP3SpecificInput7T_out src/test/resources/multiplex_test/primer_task_specific_output > test_results/diff_results/testP3SpecificInput7T_diff

echo "  -testP3SpecificInput3"
java -jar target/primer3-1.0-SNAPSHOT.jar src/test/resources/multiplex_test/primer_task_specific_input_3T -o test_results/testP3SpecificInput3_out -fasta src/test/resources/multiplex_test/targets_new_fake2.fasta

diff test_results/testP3SpecificInput3_out src/test/resources/multiplex_test/primer_task_specific_output > test_results/diff_results/testP3SpecificInput3_diff

echo "  -testP3SpecificAlbInput"
java -jar target/primer3-1.0-SNAPSHOT.jar src/test/resources/multiplex_test/primer_task_specific_alb_input -o test_results/testP3SpecificAlbInput_out -fasta src/test/resources/multiplex_test/primer_task_specific.fasta

diff test_results/testP3SpecificAlbInput_out src/test/resources/multiplex_test/primer_task_specific_output > test_results/diff_results/testP3SpecificAlbInput_diff

echo "end block TestMultiplexInput"
echo "start block TestPrimer3Main"

echo "  -testP3Case1"
echo "Ignored: Long test"
# java -jar target/primer3-1.0-SNAPSHOT.jar src/test/resources/primer1_input -o test_results/testP3Case1_out

# diff test_results/testP3Case1_out src/test/resources/primer1_output > test_results/diff_results/testP3Case1_diff

echo "  -testP3CheckInput"
java -jar target/primer3-1.0-SNAPSHOT.jar src/test/resources/primer_check_input -o test_results/testP3CheckInput_out

diff test_results/testP3CheckInput_out src/test/resources/primer_check_output > test_results/diff_results/testP3CheckInput_diff

echo "  -testP3EndPathology"
java -jar target/primer3-1.0-SNAPSHOT.jar src/test/resources/primer_end_pathology_input -o test_results/testP3EndPathology_out

diff test_results/testP3EndPathology_out src/test/resources/primer_end_pathology_output > test_results/diff_results/testP3EndPathology_diff

echo "  -testP3GCEnd"
java -jar target/primer3-1.0-SNAPSHOT.jar src/test/resources/primer_gc_end_input -o test_results/testP3GCEnd_out

diff test_results/testP3GCEnd_out src/test/resources/primer_gc_end_output > test_results/diff_results/testP3GCEnd_diff

echo "  -testP3PrimerFirstBaseIndex"
java -jar target/primer3-1.0-SNAPSHOT.jar src/test/resources/primer_first_base_index_input -o test_results/testP3PrimerFirstBaseIndex_out

diff test_results/testP3PrimerFirstBaseIndex_out src/test/resources/primer_first_base_index_output > test_results/diff_results/testP3PrimerFirstBaseIndex_diff

echo "  -testP3HighGCLoadSet"
java -jar target/primer3-1.0-SNAPSHOT.jar src/test/resources/primer_high_gc_load_set_input -o test_results/testP3HighGCLoadSet_out --p3_settings_file  src/test/resources/primer_high_tm_load_set.set

diff test_results/testP3HighGCLoadSet_out src/test/resources/primer_high_gc_load_set_output > test_results/diff_results/testP3HighGCLoadSet_diff

echo "  -testP3HighTMLoadSet"
java -jar target/primer3-1.0-SNAPSHOT.jar src/test/resources/primer_high_tm_load_set_input -o test_results/testP3HighTMLoadSet_out --p3_settings_file  src/test/resources/primer_high_gc_load_set.set

diff test_results/testP3HighTMLoadSet_out src/test/resources/primer_high_tm_load_set_output > test_results/diff_results/testP3HighTMLoadSet_diff

echo "  -testP3InternalInput"
java -jar target/primer3-1.0-SNAPSHOT.jar src/test/resources/primer_internal_input -o test_results/testP3InternalInput_out

diff test_results/testP3InternalInput_out src/test/resources/primer_internal_output > test_results/diff_results/testP3InternalInput_diff

echo "  -testP3OkRegionsInput1"
java -jar target/primer3-1.0-SNAPSHOT.jar src/test/resources/primer_ok_regions_input -o test_results/testP3OkRegionsInput1_out

diff test_results/testP3OkRegionsInput1_out src/test/resources/primer_ok_regions_output > test_results/diff_results/testP3OkRegionsInput1_diff

echo "  -testP3OkRegionsInput2"
java -jar target/primer3-1.0-SNAPSHOT.jar src/test/resources/primer_ok_regions2_input -o test_results/testP3OkRegionsInput2_out

diff test_results/testP3OkRegionsInput2_out src/test/resources/primer_ok_regions2_output > test_results/diff_results/testP3OkRegionsInput2_diff

echo "  -testP3MustMatchInput"
echo "Ignored: Long test"
# java -jar target/primer3-1.0-SNAPSHOT.jar src/test/resources/primer_must_match_input -o test_results/testP3MustMatchInput_out

#diff test_results/testP3MustMatchInput_out src/test/resources/primer_must_match_output > test_results/diff_results/testP3MustMatchInput_diff

echo "  -testP3MustOverlapPointInput"
java -jar target/primer3-1.0-SNAPSHOT.jar src/test/resources/primer_must_overlap_point_input -o test_results/testP3MustOverlapPointInput_out

diff test_results/testP3MustOverlapPointInput_out src/test/resources/primer_must_overlap_point_output > test_results/diff_results/testP3MustOverlapPointInput_diff

echo "  -testP3NumBestInput"
java -jar target/primer3-1.0-SNAPSHOT.jar src/test/resources/primer_num_best_input -o test_results/testP3NumBestInput_out

diff test_results/testP3NumBestInput_out src/test/resources/primer_num_best_output > test_results/diff_results/testP3NumBestInput_diff

echo "  -testP3OverlapJunctionInput"
java -jar target/primer3-1.0-SNAPSHOT.jar src/test/resources/primer_overlap_junction_input -o test_results/testP3OverlapJunctionInput_out

diff test_results/testP3OverlapJunctionInput_out src/test/resources/primer_overlap_junction_output > test_results/diff_results/testP3OverlapJunctionInput_diff

echo "  -testP3QualityBoundaryInput"
java -jar target/primer3-1.0-SNAPSHOT.jar src/test/resources/primer_quality_boundary_input -o test_results/testP3QualityBoundaryInput_out

diff test_results/testP3QualityBoundaryInput_out src/test/resources/primer_quality_boundary_output > test_results/diff_results/testP3QualityBoundaryInput_diff

echo "  -testP3RatInput"
java -jar target/primer3-1.0-SNAPSHOT.jar src/test/resources/primer_rat_input -o test_results/testP3RatInput_out

diff test_results/testP3RatInput_out src/test/resources/primer_rat_output > test_results/diff_results/testP3RatInput_diff

echo "  -testP3StartCodonInput"
java -jar target/primer3-1.0-SNAPSHOT.jar src/test/resources/primer_start_codon_input -o test_results/testP3StartCodonInput_out

diff test_results/testP3StartCodonInput_out src/test/resources/primer_start_codon_output > test_results/diff_results/testP3StartCodonInput_diff

echo "  -testP3HumanInput"
echo "Ignored: Very long test"
# java -jar target/primer3-1.0-SNAPSHOT.jar src/test/resources/primer_human_input -o test_results/testP3HumanInput_out

# diff test_results/testP3HumanInput_out src/test/resources/primer_human_output > test_results/diff_results/testP3HumanInput_diff

echo "  -testP3NewTaskInput"
java -jar target/primer3-1.0-SNAPSHOT.jar src/test/resources/primer_new_tasks_input -o test_results/testP3NewTaskInput_out

diff test_results/testP3NewTaskInput_out src/test/resources/primer_new_tasks_output > test_results/diff_results/testP3NewTaskInput_diff

echo "end block TestPrimer3Main"

grep -c -- "---" resourcetest_results/diff_results/* > test_results/diff_results/number_of_diff