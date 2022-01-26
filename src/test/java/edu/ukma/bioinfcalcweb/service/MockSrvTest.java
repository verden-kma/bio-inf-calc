package edu.ukma.bioinfcalcweb.service;

import edu.ukma.bioinfcalcweb.util.MockRandomProvider;
import org.junit.jupiter.api.AfterEach;
import org.junit.jupiter.api.Test;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.boot.test.context.SpringBootTest;
import org.springframework.test.context.ActiveProfiles;

import java.lang.reflect.Field;
import java.util.List;

import static org.assertj.core.api.Assertions.assertThat;

@SpringBootTest
@ActiveProfiles("test")
public class MockSrvTest {
    @Autowired
    private MockRandomProvider randomProvider;
    @Autowired
    private MotifSrv motifSrv;

    @AfterEach
    public void resetRandom() throws NoSuchFieldException, IllegalAccessException {
        Field valueIndexField = MockRandomProvider.class.getDeclaredField("randIndex");
        valueIndexField.setAccessible(true);
        valueIndexField.set(randomProvider, 0);
    }

    @Test
    public void randomMotifSearchTest() {
        String[] expected = new String[]{"TCAGTAAA", "GCGAGGTA", "AAGAAGTA", "AGGTGCAC", "ACGTGCAA"};
        List<String> actual = motifSrv.randomizedMotifSearch(new String[]{
                "CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA",
                "GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG",
                "TAGTACCGAGACCGAAAGAAGTATACAGGCGT",
                "TAGATCAAGTTTCAGGTGCACGTCGGTGAACC",
                "AATCCACCAGCTCCACGTGCAATGTTGGCCTA"
        }, 8);
        assertThat(expected).containsAll(actual);
    }

    @Test
    public void gibbsSamplerSearchTest() {
        String[] expected = new String[]{"GTGTTCAG", "GCGAGGTA", "AAGAAGTA", "AGGTGCAC", "ACGTGCAA"};
        List<String> actual = motifSrv.gibbsSamplerSearch(new String[]{
                "CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA",
                "GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG",
                "TAGTACCGAGACCGAAAGAAGTATACAGGCGT",
                "TAGATCAAGTTTCAGGTGCACGTCGGTGAACC",
                "AATCCACCAGCTCCACGTGCAATGTTGGCCTA"}, 8, 100);
        assertThat(expected).containsAll(actual);
    }
}
