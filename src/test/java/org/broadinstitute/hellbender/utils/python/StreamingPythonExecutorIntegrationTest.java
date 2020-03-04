package org.broadinstitute.hellbender.utils.python;

import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.SkipException;
import org.testng.annotations.Test;

import java.util.Collections;
import java.util.HashMap;
import java.util.Map;

public class StreamingPythonExecutorIntegrationTest extends GATKBaseTest {
    private static final String NL = System.lineSeparator();

    // NOTE: these must be kept in sync with the versions in gatkcondaenv.yml.template
    private static final Map<String, String> expectedDependencyVersions = new HashMap<String, String>() {
        private static final long serialVersionUID = 1L;
        {
            put("keras", "2.2.0");
            put("matplotlib", "2.1.0");
            put("numpy", "1.18.1");
            put("pandas", "0.21.0");
            put("scipy", "1.0.0");
            put("tensorflow", "1.12.0");
            put("theano", "1.0.4");
        }
    };

    @Test(groups = "python")
    public void testValidateDependencyVersions() {
        final StreamingPythonScriptExecutor<String> streamingPythonExecutor =
                new StreamingPythonScriptExecutor<>(true);
        Assert.assertNotNull(streamingPythonExecutor);
        Assert.assertTrue(streamingPythonExecutor.start(Collections.emptyList(), true, null));

        // validate for some key dependencies that the conda environment got the version that we asked for
        try {
            expectedDependencyVersions.forEach((p, v) -> {
                streamingPythonExecutor.sendSynchronousCommand(String.format("import %s" + NL, p));
                streamingPythonExecutor.sendSynchronousCommand(String.format("assert(%s.__version__ == '%s')" + NL, p, v));
            });
        } finally {
            streamingPythonExecutor.terminate();
        }
    }

    // This test validates that the StreamingPythonScriptExecutor will throw if the Python environment
    // is not activated, and will not pass if the environment is activated. It has to be an integration
    // test because unit tests are run on the Docker image, which has the Python environment activated.

    @Test()
    public void testRequirePythonEnvironment() {
        // This test is deliberately left out of the "python" test group in order to ensure that
        // it only executes when the Python environment has *NOT* been properly established. Also,
        // skip this test if we're running on the Docker because the Python environment is always
        // activated there.
        if (isGATKDockerContainer()) {
            throw new SkipException("Python environment validation test must be skipped when running on the Docker");
        }

        // validate that we throw if the GATK Python environment is not active
        final RuntimeException rte = Assert.expectThrows(RuntimeException.class, ()-> {
            final StreamingPythonScriptExecutor<String> streamingPythonExecutor =
                    new StreamingPythonScriptExecutor<>(PythonScriptExecutor.PythonExecutableName.PYTHON3, true);
            streamingPythonExecutor.start(Collections.emptyList());
        });

        // make sure that the underlying cause is actually a PythonScriptExecutorException
        Assert.assertEquals(rte.getCause().getClass(), PythonScriptExecutorException.class);
    }

}
