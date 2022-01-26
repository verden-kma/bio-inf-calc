package edu.ukma.bioinfcalcweb.util;

import org.springframework.context.annotation.Profile;
import org.springframework.stereotype.Component;

import java.util.Random;

@Component
@Profile("!test")
public class RealRandomProvider implements IRandomProvider {
    private final Random random = new Random();

    @Override
    public int nextInt(int upperBound) {
        return random.nextInt(upperBound);
    }

    @Override
    public double nextDouble() {
        return random.nextDouble();
    }
}
