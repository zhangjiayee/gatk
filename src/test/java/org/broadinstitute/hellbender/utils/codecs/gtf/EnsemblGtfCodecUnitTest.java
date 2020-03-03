package org.broadinstitute.hellbender.utils.codecs.gtf;

import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Test class for the ENSEMBL GTF Reader.
 * Modeled after the TableCodecUnitTest, with extras specific to this file format.
 * Created by jonn on 2020 02 28.
 */
public class EnsemblGtfCodecUnitTest extends GATKBaseTest {

    private static final String testResourceDir = publicTestDir + "org/broadinstitute/hellbender/utils/codecs/gtf/";
    private static final String eColiTestDir = publicTestDir + "org/broadinstitute/hellbender/tools/funcotator/ecoli_ds/gencode/ASM584v2/";

    @DataProvider
    public Object[][] canDecodeProvider() {

        return new Object[][] {
                { "a.tsv"     , testResourceDir, false },                                    // Wrong File name / extension
                { "a.table.gz", testResourceDir, false },                                    // Wrong File name / extension
                { "a.bed"     , testResourceDir, false },                                    // Wrong File name / extension
                { "a.bcf"     , testResourceDir, false },                                    // Wrong File name / extension
                { "a.hapmap"  , testResourceDir, false },                                    // Wrong File name / extension
                { "a.refseq"  , testResourceDir, false },                                    // Wrong File name / extension
                { "a.beagle"  , testResourceDir, false },                                    // Wrong File name / extension
                { "a.table"   , testResourceDir, false },                                    // Wrong File name / extension

                { "gencode.v26.annotation.gtf.tsv", testResourceDir, false},                 // Wrong File name / extension
                { "gencode.v26.annotation.tgz"    , testResourceDir, false},                 // Wrong File name / extension
                { "gencode.v26.annotation.tar.gz" , testResourceDir, false},                 // Wrong File name / extension

                { "gencode.gtf"                                , testResourceDir, false},    // File does not exist
                { "gencode.v26.primary_assembly.annotation.gtf", testResourceDir, false},    // File does not exist
                { "gencode.v26.long_noncoding_RNAs.gtf"        , testResourceDir, false},    // File does not exist

                { "gencode.invalid_short_header.gtf"           , testResourceDir, false},    // File exists, has invalid header
                { "gencode.invalid_malformed_header.gtf"       , testResourceDir, false},    // File exists, has invalid header
                { "gencode.invalid_malformed_header_desc.gtf"  , testResourceDir, false},    // File exists, has invalid header
                { "gencode.invalid_malformed_header_prov.gtf"  , testResourceDir, false},    // File exists, has invalid header
                { "gencode.invalid_malformed_header_cont.gtf"  , testResourceDir, false},    // File exists, has invalid header
                { "gencode.invalid_malformed_header_form.gtf"  , testResourceDir, false},    // File exists, has invalid header
                { "gencode.invalid_malformed_header_date.gtf"  , testResourceDir, false},    // File exists, has invalid header

                { "gencode.valid1.gtf"                           , testResourceDir, false},   // Not an Ensembl GTF file
                { "gencode.valid_gencode_file2.gtf"              , testResourceDir, false},   // Not an Ensembl GTF file
                { "gencode.and.this.is.a.valid.one.too.table.gtf", testResourceDir, false},   // Not an Ensembl GTF file

                { "Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.44.gtf", eColiTestDir, true},   // Name doesn't start with 'gencode'
        };
    }

    @Test(dataProvider = "canDecodeProvider")
    public void testCanDecode(final String fileName, final String containingFolder, final boolean expected) {
        final EnsemblGtfCodec ensemblGtfCodec = new EnsemblGtfCodec();
        Assert.assertEquals(ensemblGtfCodec.canDecode(containingFolder + fileName), expected, fileName);
    }

    @DataProvider
    public Object[][] headerProvider() {
        return new Object[][] {

                { new ArrayList<String>(), false },                             // Wrong length header
                { Arrays.asList( "",
                        "",
                        "",
                        "",
                        ""  ),
                        false },                                // Bad content
                { Arrays.asList( "##descr",
                        "##provider: GENCODE",
                        "##contact: gencode-help@sanger.ac.uk",
                        "##format: gtf",
                        "##date: 2017-04-08" ),
                        false },                                // Bad header - description
                { Arrays.asList( "##description: THIS IS A SAMPLE",
                        "##provider: GARBAGEDAY",
                        "##contact: gencode-help@sanger.ac.uk",
                        "##format: gtf",
                        "##date: 2017-04-08" ),
                        false },                                // Bad header - provider
                { Arrays.asList( "##description: THIS IS A SAMPLE",
                        "##provider: GENCODE",
                        "##contact: gencode@NORTHPOLE.pl",
                        "##format: gtf",
                        "##date: 2017-04-08" ),
                        false },                                // Bad header - contact
                { Arrays.asList( "##description: THIS IS A SAMPLE",
                        "##provider: GENCODE",
                        "##contact: SANTACLAUSE@sanger.ac.uk",
                        "##format: gtf",
                        "##date: 2017-04-08" ),
                        false },                                // Bad header - contact
                { Arrays.asList( "##description: THIS IS A SAMPLE",
                        "##provider: GENCODE",
                        "##contact: gencode-help@sanger.ac.uk",
                        "##format: dumpy",
                        "##date: 2017-04-08" ),
                        false },                                // Bad header - format
                { Arrays.asList( "##description: THIS IS A SAMPLE",
                        "##provider: GENCODE",
                        "##contact: gencode-help@sanger.ac.uk",
                        "##format: gtf",
                        "##doom: ID Software" ),
                        false },                                // Bad header - date
                { Arrays.asList( "##description: evidence-based annotation of the human genome (GRCh37), version 19 (Ensembl 74)",
                        "##provider: GENCODE",
                        "##contact: gencode@sanger.ac.uk",
                        "##format: gtf",
                        "##date: 2014-07-25" ),
                        false },                                // Good GENCODE Header, but BAD ENSEMBL header!
                { Arrays.asList( "##description: evidence-based annotation of the human genome (GRCh38), version 26 (Ensembl 88)",
                        "##provider: GENCODE",
                        "##contact: gencode-help@sanger.ac.uk",
                        "##format: gtf",
                        "##date: 2014-07-25" ),
                        false },                                 // Good GENCODE Header, but BAD ENSEMBL header!

                // -------------

                { Arrays.asList( "#!genome-build ASM584v2",
                        "#!genome-version ASM584v2",
                        "#!genome-date 2014-08",
                        "#!genome-build-accession GCA_000005845.2",
                        "#!genebuild-last-updated 2014-08" ),
                        true },                                           // Good Ensembl GTF Header!

                { Arrays.asList( "ASM584v2",
                        "#!genome-version ASM584v2",
                        "#!genome-date 2014-08",
                        "#!genome-build-accession GCA_000005845.2",
                        "#!genebuild-last-updated 2014-08" ),
                        false },                                           // Bad header - genome-build
                { Arrays.asList( "#!genome-build ASM584v2",
                        "ASM584v2",
                        "#!genome-date 2014-08",
                        "#!genome-build-accession GCA_000005845.2",
                        "#!genebuild-last-updated 2014-08" ),
                        false },                                           // Bad header - genome-version
                { Arrays.asList( "#!genome-build ASM584v2",
                        "#!genome-version ASM584v2",
                        "#2014-08",
                        "#!genome-build-accession GCA_000005845.2",
                        "#!genebuild-last-updated 2014-08" ),
                        false },                                           // Bad header - genome-date
                { Arrays.asList( "#!genome-build ASM584v2",
                        "#!genome-version ASM584v2",
                        "#!genome-date 2014-08",
                        "#GCA_000005845.2",
                        "#!genebuild-last-updated 2014-08" ),
                        false },                                           // Bad header - genome-build-accession
                { Arrays.asList( "#!genome-build ASM584v2",
                        "#!genome-version ASM584v2",
                        "#!genome-date 2014-08",
                        "#!genome-build-accession GCA_000005845.2",
                        "#2014-08" ),
                        false },                                           // Bad header - genebuild-last-updated

        };
    }

    @Test(dataProvider = "headerProvider")
    public void testValidateHeader(final List<String> header, final boolean expected ) {
        final EnsemblGtfCodec ensemblGtfCodec = new EnsemblGtfCodec();
        Assert.assertEquals( ensemblGtfCodec.validateHeader(header), expected );
    }

}
