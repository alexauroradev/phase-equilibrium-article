#include "PhaseSolver.h"
#include <iostream>

int main() {
	la::vector<double> c(3);

	c[0] = 0.5; /* H2O */
	c[1] = 0.25; /* C5H12 */
	c[2] = 0.25; /* C10H22 */
	//c[3] = 0.01; /* C15H32 */
	//c[4] = 0.08; /* N2 */
	//c[5] = 0.005; /* CO2 */
	//c[6] = 0.02; /* O2 */
	//c[7] = 1; /* Ske */
	//c[8] = 0.01; /* Ker */

	double p = 1e7;
	double Tapprox = 300;

	PhaseSolver ps(/* LO */ 2,/* HO */ 0,/* G */ 0,/* S */ 0, /* eps */ 1e-6);
	//PhaseSolver ps(/* LO */ 2,/* HO */ 1,/* G */ 3,/* S */ 2, /* eps */ 1e-6);
	ps.setEquilibriumConstantLog(0, new Kappa(11.8572e9, 3816.44, 46.13)); /* H2O */
	ps.setEquilibriumConstantLog(1, new Kappa(1.0026e9, 2477.07, 39.94)); /* C5H12 */
	ps.setEquilibriumConstantLog(2, new Kappa(1.1981e9, 3456.80, 78.67)); /* C10H22 */

	ps.setStandartEnthalpy(0, new GasEnthalpy(31.876e3, 6.493)); /* H2O */
	ps.setStandartEnthalpy(1, new GasEnthalpy(6.292e3, 367.4)); /* C5H12 */
	ps.setStandartEnthalpy(2, new GasEnthalpy(12.291, .7160)); /* C10H22 */

	//ps.setStandartEnthalpy(3, new GasEnthalpy(2e3, 0)); /* C15H32 */

	//ps.setStandartEnthalpy(4, new GasEnthalpy(30.288e3, -2.572)); /* N2 */
	//ps.setStandartEnthalpy(5, new GasEnthalpy(16.864e3, 106.3)); /* CO2 */
	//ps.setStandartEnthalpy(6, new GasEnthalpy(27.627e3, 6.437)); /* O2 */

	//ps.setStandartEnthalpy(7, new GasEnthalpy(209.6e3, 0)); /* Ske, 800 J/kg-K * 262 kg/kmol */
	//ps.setStandartEnthalpy(8, new GasEnthalpy(1.5e6, 0)); /* Ker, 1500 J/kg-K * 1000 kg/kmol */

	ps.setVaporizationEnthalpy(0, new VaporizationEnthalpy(4.82e6, 647.3));
	ps.setVaporizationEnthalpy(1, new VaporizationEnthalpy(3.745e6, 469.6));
	ps.setVaporizationEnthalpy(2, new VaporizationEnthalpy(5.579e6, 617.6));

	ps.setTempRange(100, 2000); /* Range for which correlations are valid */

	ps.setParams(c, p, Tapprox, 1e7);
	std::cerr << "Test F : \n" << (ps.checkF() ? "OK" : "Failed") << std::endl;
	std::cerr << "Test Fx: \n" << (ps.checkFx() ? "OK" : "Failed") << std::endl;
	std::cerr << "Test Fa: \n" << (ps.checkFa() ? "OK" : "Failed") << std::endl;

	std::cout << "eta,ome,lam0,lam1,T" << std::endl;
	//std::cout << "eta,ome,lam0,lam1,T,dode,dl0de,dl1de,dtde" << std::endl;

	for (double eta = 1e5; eta < 3.5e7; eta += 1000) {
		ps.setParams(c, p, Tapprox, eta);

		la::vector<double> center = ps.solve();
		la::matrix<double> der = ps.getDerivatives();

		double res;

		Status s = ps.getStatus(res);

		if (s != Ok)
			std::cerr << s << " " << res << std::endl;

		std::cout << util::format("%e,%.20e,%.20e,%.20e,%.20e\n",
			eta,
			center[0], center[1], center[2], center[3]//,
			//der(0, 1), der(1, 1), der(2, 1), der(3, 1)
			);
	}
	return 0;
}
