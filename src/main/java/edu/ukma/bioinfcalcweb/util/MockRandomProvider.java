package edu.ukma.bioinfcalcweb.util;

import org.springframework.context.annotation.Profile;
import org.springframework.stereotype.Component;

import java.util.Arrays;
import java.util.List;

@Component
@Profile("test")
public class MockRandomProvider implements IRandomProvider {
    private static final List<Double> RAND_VALS = Arrays.asList(
            0.730967787376657,
            0.24053641567148587,
            0.6374174253501083,
            0.5504370051176339,
            0.5975452777972018,
            0.3332183994766498,
            0.3851891847407185,
            0.984841540199809,
            0.8791825178724801,
            0.9412491794821144,
            0.27495396603548483,
            0.12889715087377673,
            0.14660165764651822,
            0.023238122483889456,
            0.5467397571984656,
            0.9644868606768501,
            0.10449068625097169,
            0.6251463634655593,
            0.4107961954910617,
            0.7763122912749325);

    private int randIndex = 0;

    @Override
    public int nextInt(int top) {
        double v = RAND_VALS.get(randIndex++) * top;
        randIndex %= RAND_VALS.size();
        return (int) v;
    }

    @Override
    public double nextDouble() {
        double v = RAND_VALS.get(randIndex++);
        randIndex %= RAND_VALS.size();
        return v;
    }
}
