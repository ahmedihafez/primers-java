Primer3
========================
# Building
For build primer3 with maven run:

`mvn clean package`

make sure that you are in the same folder than `pom.xml`.

# Testing

Use the `test.sh` to run test the builded jar.

The outputs are stored in `test_results`, each output is compared to a reference and the diferences are stored in `test_results/diff_results`.
