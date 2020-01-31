package org.broadinstitute.hellbender.tools;

import htsjdk.samtools.SamReaderFactory;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.SamAssertionUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;

public final class PrintReadsIntegrationTest extends AbstractPrintReadsIntegrationTest {

    @Test(dataProvider="testingData")
    public void testFileToFileWithMD5(String fileIn, String extOut, String reference) throws Exception {
        doFileToFile(fileIn, extOut, reference, true);
    }

    @Test
    public void testNoConflictPG() throws IOException {
        final File inFile = new File(TEST_DATA_DIR, "print_reads_withPG.sam");
        final File outFile = GATKBaseTest.createTempFile("testNoConflictRG", ".sam");
        final String[] args = new String[] {
                "--input" , inFile.getAbsolutePath(),
                "--" + StandardArgumentDefinitions.ADD_OUTPUT_SAM_PROGRAM_RECORD,
                "--output", outFile.getAbsolutePath()
        };
        runCommandLine(args);

        //Make sure contents are the same
        SamAssertionUtils.assertSamsEqual(outFile, inFile);

        //input has GATK PrintReads not NOT GATK PrintReads.1 in headers
        Assert.assertNotNull(SamReaderFactory.makeDefault().open(inFile).getFileHeader().getProgramRecord("GATK PrintReads"));
        Assert.assertNull(SamReaderFactory.makeDefault().open(inFile).getFileHeader().getProgramRecord("GATK PrintReads.1"));

        //output has both GATK PrintReads and GATK PrintReads.1 in headers
        Assert.assertNotNull(SamReaderFactory.makeDefault().open(outFile).getFileHeader().getProgramRecord("GATK PrintReads"));
        Assert.assertNotNull(SamReaderFactory.makeDefault().open(outFile).getFileHeader().getProgramRecord("GATK PrintReads.1"));
    }

    @Test
    public void testHTTP(){
        ArgumentsBuilder args = new ArgumentsBuilder();
        File out = createTempFile("out", ".bam");
        args.addArgument("I", "https://storage.googleapis.com/hellbender/test/resources/benchmark/CEUTrio.HiSeq.WEx.b37.NA12892.bam")
                        .addArgument("read-index", "https://storage.googleapis.com/hellbender/test/resources/benchmark/CEUTrio.HiSeq.WEx.b37.NA12892.bam.bai")
                        .addOutput(out)
                        .addInterval(new SimpleInterval("3",1_000_000, 1_000_001))
                        .addInterval(new SimpleInterval("3", 1_000_003, 1_000_100))
                        .addInterval(new SimpleInterval("20", 1_000_000, 1_100_000));

        runCommandLine(args);


    }
}