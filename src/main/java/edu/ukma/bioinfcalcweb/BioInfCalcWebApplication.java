package edu.ukma.bioinfcalcweb;

import edu.ukma.bioinfcalcweb.controller.IControllerPackageMarker;
import edu.ukma.bioinfcalcweb.service.IServicePackageMarker;
import edu.ukma.bioinfcalcweb.util.UtilsPackageMarker;
import org.springframework.boot.SpringApplication;
import org.springframework.boot.autoconfigure.SpringBootApplication;
import org.springframework.context.annotation.ComponentScan;

@SpringBootApplication
@ComponentScan(basePackageClasses = {IControllerPackageMarker.class, IServicePackageMarker.class, UtilsPackageMarker.class})
public class BioInfCalcWebApplication {
    // "@RestController" is not user because of Tro
    public static void main(String[] args) {
        SpringApplication.run(BioInfCalcWebApplication.class, args);
    }
}
