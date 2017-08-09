
#include <math.h>
#include "IsentropFlow.h"
#include "SIE_Math.h"

double VATfM(const double gamma, const double Mach) {
        return Mach / sqrt(1 + (gamma - 1) / 2 * pow(Mach,2));
}

double VELfTTM(const double gamma, const double Tt, const double Mach) {
		double Ts = Tt / TRfM(gamma, Mach);
		double VsonicRT = sqrt(gamma * Rgas * Ts);
        return VATfM(gamma, Mach) * VsonicRT;
}

double QoverPT(const double gamma, const double Mach) {
        //   Computes Ratio of Dynamic Pressure to Total Pressure given Gamma and Mach Number
        //       Note: Dynamic Pressure, Q = 0.5*RHO*V^2
        //   gamma   ~Specific Heat Ratio
        return 0.5 * gamma * pow(Mach,2) * (1 / PRfM(gamma, Mach));
}

double PMturn_NUfM(const double gamma, const double Mach) {
        //Calculate Prandtl-Meyer turn angle, nu, as function of Mach
        double Nu;
		double gp1, gm1;
        gp1 = gamma + 1;
        gm1 = gamma - 1;
        Nu = sqrt(gp1 / gm1) * atan(sqrt(gm1 / gp1 * (pow(Mach,2) - 1))) - atan(sqrt(pow(Mach,2) - 1));
        return Nu * 180 / PI;
}

double PMturn_MfNU(const double gamma, const double Nu) {
        //Calculate Prandtl-Meyer expansion Mach as a function of nu

        double M1, dM, kdir, tol, Nu1;

        M1 = 1.0001;
        dM = 0.05;
        kdir = 1;
        tol = 0.00001;

        int i;
        for (i = 1; i <= 100; i++) {
            Nu1 = PMturn_NUfM(gamma, M1);
            if (Nu1 - Nu > tol) {  //M too big
                if (kdir == 1) {
                    kdir = -1;
                    dM = dM / 2;
                }
            } else if (Nu - Nu1 > tol) { //M too small
                if (kdir == -1) {
                    kdir = 1;
                    dM = dM / 2;
                }
            } else {
                break;
            }
            M1 = M1 + kdir * dM;
        }

        return M1;

}

double MfNu(const double gamma, const double Nu) {

        // Determines Mach number from Prandtl-Meyer turn angle.
        // Accurate to atleast turn angles of 95 degrees
        // written by DAH

        const double Macc = 0.00001;
        const int jmax = 1000;

        double M1, M2, M3, dM, dNu1, dNu2, dNu3;
        M1 = 1;
        M2 = 100;

        dNu1 = PMturn_NUfM(gamma, M1) - Nu;
        dNu2 = PMturn_NUfM(gamma, M2) - Nu;

        int j = 0;

        while (j <= jmax) {
            M3 = 0.5 * (M1 + M2);
            dNu3 = PMturn_NUfM(gamma, M3) - Nu;
            if (dNu3 * dNu1 <= 0) {
                dM = M3 - M2;
                M2 = M3;
                dNu2 = dNu3;
            } else {
                dM = M3 - M1;
                M1 = M3;
                dNu1 = dNu3;
            }
            if ((fabs(dM) < Macc) || (dNu2 == 0)) {
            	j = jmax;
            	j = j + 1;
            }
        }

        return M2;

}

double Mu1(const double Mach) {
        //Mach angle relative to Prandtl-Meyer Turn
        return 90 - (acos(1 / Mach)) * 180 / PI;
}

double Mu2(const double gamma, const double Mach) {
        //Mach angle relative to wall upstream of Prandtl-Meyer Turn
        //Positive angles are into the approach flow
        //Negative angles are away from the approach flow

        double Mux, Nux;
        Mux = Mu1(Mach);
        
        double gp1, gm1;
        gp1 = gamma + 1;
        gm1 = gamma - 1;
        Nux = sqrt(gp1 / gm1) * atan(sqrt(gm1 / gp1 * (pow(Mach,2) - 1))) - atan(sqrt(pow(Mach,2) - 1));
        return Mux - (Nux * 180 / PI);
}

double P2P1_PMexpansion(const double gamma, const double Mach, const double delta_ang_deg) {
        //   gamma= ratio of specific heats
        //   Mach = approach Mach number
        //   nu   = incremental change in angle (degrees)

        double Nu, PsPtx, dNu, term1, term2, term3, term4;
        Nu = delta_ang_deg * PI / 180;
        PsPtx = 1;
        dNu = Nu / 100;

        int i;
        for (i = 1; i <= 100; i++) {
            term1 = 1 - (gamma * pow(Mach,2)) / sqrt((pow(Mach,2) - 1)) * dNu;
            term2 = gamma * pow(Mach,2) * ((gamma + 1) * pow(Mach,4) - 4 * (pow(Mach,2) - 1)) / (4 * pow(pow(Mach,2) - 1, 2) * pow(dNu,2));
            term3 = gamma * pow(Mach,2) / (2 * pow((pow(Mach,2) - 1),3.5)) * pow(dNu,3);
            term4 = (((gamma + 1) / 6 * pow(Mach,8)) 
                    - ((5 + 7 * gamma - 2 * pow(gamma,2)) / 6 * pow(Mach,6)) 
                    + (5 / 3 * (gamma + 1) * pow(Mach,4) - 2 * pow(Mach,2) + 4 / 3));
            PsPtx = (term1 + term2 - term3 * term4) * PsPtx;
        }
        return PsPtx;
}

/*
// JGB: such an easy equation for sound speed in SI units, this equation is commented out
double VsonicRT(const double gamma, const double Ts, const double Rgas) {
        // Computes Speed of Sound in a Gas Given Gamma, Gas Constant, and Absolute Temperature
        //   gamma       ~specific heat ratio for gas
        //   Rgas        ~gas constant (ft^2/(s^2*R))
        //   Ts          ~Absolute Temperature (deg. R)
        //   VsonicRT    ~Sonic Velocity (ft/sec)

        if (Constants.English.Grav * gamma * Rgas * Ts < 0) {
            return -1000000000000.0;
        } else {
            return sqrt(Constants.English.Grav * gamma * Rgas * Ts)
        }

}
*/

double TRfM(const double gamma, const double Mach) {
        //   Computes Total-to-Static Temperature ratio (Ttotal/Tstatic) for Isentropic Flow (given Gamma and Mach Number)
        //   gamma   ~Specific Heat Ratio
        return 1 + ((gamma - 1) / 2) * pow(Mach,2);
}

double MfTR(const double gamma, const double Tr) {
        //   Computes Isentropic Mach Number given Gamma and Total-to-Static Temperature Ratio
        //   gamma   ~Specific Heat Ratio
        return sqrt((2 / (gamma - 1)) * (Tr - 1));
}

double PRfM(const double gamma, const double Mach) {
        //   Computes Total to Static Pressure ratio (Ptotal/Pstatic) for Isentropic Flow (given Gamma and Mach Number)
        //   gamma   ~Specific Heat Ratio
        return pow((1 + ((gamma - 1) / 2) * pow(Mach,2)), (gamma / (gamma - 1)));
}

double MfPR(const double gamma, const double Pr) {
        //   Computes Isentropic Mach Number given Gamma and Total-to-Static Pressure Ratio
        //   gamma   ~Specific Heat Ratio
        return sqrt(2 / (gamma - 1) * (pow(Pr, (gamma - 1) / gamma) - 1));
}

double DRfM(const double gamma, const double Mach) {
        //   Computes Total to Static Density ratio for Isentropic Flow (given Gamma and Mach Number)
        //   gamma   ~Specific Heat Ratio
        return pow((1 + ((gamma - 1) / 2) * pow(Mach,2)), (1 / (gamma - 1)));
}

double MfDR(const double gamma, const double Dr) {
        //   Computes Isentropic Mach Number given Gamma and Total-to-Static Density Ratio
        //   gamma   ~Specific Heat Ratio
        return sqrt((2 / (gamma - 1)) * (pow(Dr,(gamma - 1)) - 1));
}

double FPfM(const double gamma, const double Mach) {
        //  Computes Isentropic Flow Parameter (W*sqrt(Ttotal)/(Ptotal*A)) given Gamma and Mach Number
        double a = 1 + (gamma-1)/2 * pow(Mach,2);
        double b = pow(a, (gamma + 1)/(2*(1-gamma)));
        double fp = Mach * sqrt(gamma/Rgas) * b;
		return fp;
}

double MfFP(const double gamma, const double wtap, const FlowRegime regime) {
        //  THIS ROUTINE ITERATES ON MACH NUMBER AND CALCULATES PRESSURE RATIO
        //  GIVEN A FLOW PARAMETER

        double tol = 0.000001;
        int iMax = 200;
        double Rmach, RM;

        double gp1 = gamma + 1;
        double gm1 = gamma - 1;

        if (regime == kFlowSubsonic) RM = 0.1;
        if (regime == kFlowSupersonic) RM = 2.0;

        Rmach = RM;
        int i;
        for (i=1; i<=iMax; i++) {
            if (regime == kFlowSubsonic) {
                Rmach = wtap * pow((1.0 + gm1 / 2.0 * pow(RM, 2)), (gp1 / 2.0 / gm1));
            } else {
                Rmach = sqrt((pow((wtap / Rmach), (2 * (1 - gamma) / gp1) - 1)) * 2 / gm1);
            }
            if (fabs((RM - Rmach) / RM) <= tol) { 
            	break;
            } else {
            	RM = Rmach;
            }
        }
        return RM;
}

double ARfM(const double gamma, const double Mach) {
        //   Computes Ratio of Area to Sonic Area (A/A*) for Isentropic Flow (given Gamma and Mach Number)
        //   gamma   ~Specific Heat Ratio
        double gm1 = gamma - 1;
        double gp1 = gamma + 1;
        return (pow(((1 + 0.5 * gm1 * pow(Mach,2)) / (gp1 / 2)), (gp1 / (2 * gm1)))) / Mach;
}

double ARfM_brent(double Mach, ...) {
	
	//declare the variable argument list
	va_list arg_list;

	//initialize arg_list
	va_start(arg_list, Mach);

	//only get the next argument in the va_list for the ARfM function any arguments after that
	// are disregarded
	double gamma = va_arg(arg_list, double);
	double arGuess = ARfM(gamma, Mach);

	//stop any remaining calls to va_arg
	va_end(arg_list);
	
	//return the area ratio for the given Mach number
	return arGuess;

}
double MfAR(const double gamma, const double Aratio, const FlowRegime flowSpeed) {

	double M1guess, M2guess;
	if (flowSpeed == kFlowSubsonic) {
		M1guess = 0.0000001;
		M2guess = 1.0;
	} else if (flowSpeed == kFlowSupersonic) {
		M1guess = 1.0;
		M2guess = 1000;
	}

	BrentFunc fBrent = &ARfM_brent;
	const double tol = 1e-6;
	const double fBrentArgs[] = {gamma};
	double rsltMach = Brent(M1guess, M2guess, tol, fBrent, fBrentArgs);
	
	return rsltMach;
	
}

/*
double MfAR(const double gamma, const double ARatio) {
        //   Computes Mach number as a function of A/A*
        //     if (A/A*>0 returns supersonic solution
        //     if (A/A*<0 returns subsonic solution
        //
        double A[20], M[20];
        double AR, XNZ, XAZ;
        int mySign, i;
        double tol = 0.0001;

        double gp1 = gamma + 1;
        double gm1 = gamma - 1;
        double g1 = 2 / gp1;
        double g2 = gm1 / 2;
        //
        mySign = ARatio / fabs(ARatio);
        AR = fabs(ARatio);
        if (AR <= 1) {
            return 1.0;
        } else {
            M[0] = 1 + 0.5 * mySign;
            XNZ = 3;
            for (i = 1; i <= 20; i++) {
                A[i] = 1 / M[i] * (g1 * (1 + g2 * pow(M[i],2))) ^ (gp1 / (2 * gm1)));
                if (fabs(AR - A[i]) / AR <= tol) {
                    break;
                }
                if (i > 1) {
                    XNZ = log((A[i - 1] - 1) / (A[i] - 1)) / log((M[i - 1] - 1) / (M[i] - 1));
                }
                XAZ = (A[i] - 1) / pow(fabs(M[i] - 1), XNZ);
                M[i + 1] = pow((AR - 1) / XAZ),(1 / XNZ)) * mySign + 1;
                if (M[i + 1] <= 0) {
                    M[i + 1] = 0.01 / i;
                }
            }
            return M[i]
        }
}
*/

/*
double FPSfM(const double gamma, const double Mach) {
        //   Computes Isentropic Flow Parameter (W*sqrt(T)/(Pstatic*A)) given Gamma, Gas Constant, and Mach Number
        //   gamma   ~Specific Heat Ratio
        //   R ~ Gas Constant (ft^2/(s^2*R))
        return (((Constants.English.Grav * gamma / (Rgas / Constants.English.Grav)) ^ 0.5) * Mach * (1 + ((gamma - 1) / 2) * pow(Mach,2)) ^ ((gamma + 1) / (2 * (1 - gamma)))) * ((1 + ((gamma - 1) / 2) * pow(Mach,2)) ^ (gamma / (gamma - 1)))
}

double MfFPsub(ByVal gamma As Double, ByVal fp As Double, Optional ByVal Rgas As Double = 53.35) {
        //   Computes Isentropic Mach Number given Gamma, Gas Constant, and Flow Parameter (W*sqrt(T)/(Ptotal*A))
        //       Note: For Subsonic Flow only
        //   gamma   ~Specific Heat Ratio
        //   R ~ Gas Constant (ft^2/(s^2*R))
        //   fp  ~Flow Parameter ((lbm/s)*(deg.R)^0.5)/(psi*in^2)

        Dim fpchoked As Double, Mold As Double, Mnew As Double, i As Integer
        Const tol1 As Double = 0.00005
        Const tol2 As Double = 0.0000001
        Const iMax As Integer = 500
        Dim Test As Double

        if (fp = 0) {
            MfFPsub = 0
            Exit Function
        }

        Test = sqrt(fp)
        fpchoked = FPfM(gamma, 1)
        if (fabs(fp - fpchoked) < tol1 Or fp > fpchoked) { return 1

        Mold = 0.5
        i = 0
10:     REM
        Mnew = fp / (FPfM(gamma, Mold) / Mold)
        i = i + 1
        if (fabs(Mnew - Mold) < tol2) { return Mnew
        if (i > iMax) { Test = sqrt(iMax - i) //Flag to indicate too many iterations
        Mold = Mnew
        GoTo 10

}

double MfFPS(ByVal gamma As Double, ByVal fps As Double, Optional ByVal Rgas As Double = 53.35) {
        //   Computes Isentropic Mach Number given Gamma, Gas Constant, and Static Flow Parameter (W*sqrt(T)/(Pstatic*A))
        //   gamma   ~Specific Heat Ratio
        //   R ~ Gas Constant (ft^2/(s^2*R))
        //   fps  ~((lbm/s)*(deg.R)^0.5)/(psi*in^2)

        Dim Mold As Double, Mnew As Double, i As Integer
        Const tol As Double = 0.0000001
        Const iMax As Integer = 500
        Dim Test As Double

        if (fps = 0) {
            return 0
        }

        Mold = 0.5
        i = 1
10:     REM
        Mnew = fps / (FPfM(gamma, Mold) * PRfM(gamma, Mold) / Mold)
        i = i + 1
        if (fabs(Mnew - Mold) < tol) { return Mnew
        if (i > iMax) { Test = sqrt(iMax - i) //Flag to indicate too many iterations
        Mold = Mnew
        GoTo 10

}
*/

/*
double MfQPTsub(ByVal gamma As Double, ByVal QPT As Double) {
        //   Computes a Subsonic Mach Number as function of Dynamic Pressure to Total Pressure Ratio (QoverPT) and Gamma
        //   gamma   ~Specific Heat Ratio

        Dim Mold As Double, Mnew As Double, i As Integer
        Const tol As Double = 0.0000001
        Const iMax As Integer = 500
        Dim Test As Double

        if (QPT = 0) {
            return 0
        }

        Mold = 0.5
        i = 1
10:     REM
        Mnew = sqrt(2.0 * QPT * PRfM(gamma, Mold) / gamma);
        i = i + 1
        if (fabs(Mnew - Mold) < tol) { return Mnew
        if (i > iMax) { Test = sqrt(iMax - i) //Flag to indicate too many iterations
        Mold = Mnew
        GoTo 10

}
*/

double FSAOAS(const double gamma, const double Aratio) {
        //   Computes Stream Thrust Coefficient  as a function of A/A* and gamma
        //   gamma   ~Specific Heat Ratio
        //   ARATIO  ~Area Ratio, A/A*
        double Mach = MfAR(gamma, Aratio, kFlowSupersonic);
        double PsPt = 1 / PRfM(gamma, Mach);
        return PsPt * Aratio * (1 + gamma * pow(Mach,2));
}

double FOPTA(const double gamma, const double NPR) {
        //   Computes Ideal Thrust Coefficient  as a function of NPR and gamma
        //   gamma   ~Specific Heat Ratio
        //   NPR ~Nozzle Pressure
        double gp1 = gamma + 1;
        double gm1 = gamma - 1;
        return sqrt(2 * pow(gamma,2) / gm1 * pow((2 / gp1), (gp1 / gm1)) * (1 - pow((1 / NPR), (gm1 / gamma))));
}

double FOWT(const double gamma, const double NPR) {
        //   Computes Ideal Thrust Coefficient (F/(W*TT^0.5)) as a function of NPR and gamma
        double gm1 = gamma - 1;
        return sqrt(2 * gamma * Rgas / gm1 * (1 - pow((1 / NPR),(gm1 / gamma))));
}

double Cv(const double Cs, const double gamma, const double AeOAs, const double NPR) {
        //   Computes Ideal Thrust Coefficient  as a function of NPR, Cs, A/A*, and gamma
        //   gamma   ~Specific Heat Ratio
        //   NPR ~Nozzle Pressure
        return (Cs * FSAOAS(gamma, AeOAs) - AeOAs / NPR) / FOPTA(gamma, NPR);
}

double CvfCs(const double Cs, const double Cd, const double gamma
        , const double ARatio, const double NPR) {
        //   Computes Ideal Thrust Coefficient  as a function of NPR, Cs, A/A*, and gamma

        int mySign = 1;
        double ARCd = ARatio / Cd;
        if (ARatio == 1) ARCd = 1;
        //Lines added to account for unchoked nozzle ==> (fsA/A*) = f(NPR)
        double PRchoke1 = PRfM(gamma, MfAR(gamma, ARCd, kFlowSupersonic));
        //double PRchoke2 = PRfM(gamma, 1);
        double PRchoke = PRchoke1;
        if (NPR < PRchoke) {
            ARCd = ARfM(gamma, MfPR(gamma, NPR));
            mySign = -1;
        }

        return (Cs * FSAOAS(gamma, mySign * ARCd) - ARCd / NPR) / FOPTA(gamma, NPR);
}

double CsfCv(const double CV, const double Cd, const double gamma
	, const double ARatio, const double NPR) {

        //   Computes Ideal Thrust Coefficient  as a function of NPR, Cs, A/A*, and gamma

        int mySign = 1;
        double ARCd = ARatio / Cd;
        if (ARatio == 1) ARCd = 1;
        //Lines added to account for unchoked nozzle ==> (fsA/A*) = f(NPR)
        double PRchoke1 = PRfM(gamma, MfAR(gamma, ARCd, kFlowSupersonic));
        //double PRchoke2 = PRfM(gamma, 1);
        double PRchoke = PRchoke1;
        if (NPR < PRchoke) {
            ARCd = ARfM(gamma, MfPR(gamma, NPR));
            mySign = -1;
        }

        return (CV * FOPTA(gamma, NPR) + ARCd / NPR) / FSAOAS(gamma, mySign * ARCd);
}


/*
 * JGB: following routines calculate mach number from area ratio using root finding functions
 * Not implemented yet until all root finding functions are defined and tested/validated
double MachfAratio(const double gamma, const double ARatio, const FlowRegime regime) {

        double MachMin, MachMax, Mach;
        double tol = 0.000000001;
        Dim AL As New ArrayList
        AL.Add(gamma)
        AL.Add(ARatio)

        Select Case FlowRegime
            Case Regime.Subsonic
                MachMin = 0.00001
                MachMax = 1
            Case Regime.Supersonic
                MachMin = 1
                MachMax = 15 //We should change this to a secant solver because the upper bound is unknown
        End Select

        //Aerothermal.MathFuncs.PlotFunction(AddressOf CalcDeltaAR, MachMin, MachMax, 1000, CType(AL, ArrayList))
        Mach = Brent(MachMin, MachMax, tol, AddressOf CalcDeltaAR, CType(AL, ArrayList))

        return Mach

}
 double CalcDeltaAR(const double MachGuess, void data) {

        Dim gamma As Double = CType(AL(0), Double)
        Dim Aratio As Double = CType(AL(1), Double)

        Dim ARguess As Double = ARATIOfM(gamma, MachGuess)

        Dim DeltaAR As Double = Aratio - ARguess

        return DeltaAR
}
*/

/*
double VsonicGamfT(ByVal T As Double) {
        Dim gamma As Double
        // Computes Speed of Sound in Air Given Absolute Temperature
        //   Vsonic  ~Sonic Velocity (ft/sec)
        //   T       ~Absolute Temperature (deg. R)
        //   Note: Valid in range 509 to 3060 degrees R

        if (T < 509 Or T > 3060) {
            MsgBox("Error: Temperature Out of Range (Gamma Correlation)")
            Exit Function
        }

        gamma = 4.825598E-18 * T ^ 5 - 0.00000000000004997829 * T ^ 4 + 0.0000000001961552 * T ^ 3 - 0.000000352011 * T ^ 2 + 0.0002303887 * T + 1.351311
        return VsonicRT(gamma, T)
}
*/

double CD_ASME(const double ReD) {
        //flow coefficient for ASME nozzle
        //Fundamentals of Temperature, Pressure, and Flow Measurements
        //Robert Benedict, pg 429 (provided by Peter Giese (ASE))
        double Cd, Cd1, Cd2;

        if ((ReD >= 1000) && (ReD <= 200000)) { //laminar solution
            Cd = 1 - 6.92 / sqrt(ReD);
        } else if (ReD >= 1000000 && ReD <= 10000000) { //turbulent solution
            Cd = 1 - 0.184 / pow(ReD, 0.2);
        } else if ((ReD > 200000) && (ReD < 1000000)) { //REd is out of range of correlations
            Cd2 = 1 - 0.184 / pow(ReD, 0.2);
            Cd1 = 1 - 6.92 / sqrt(ReD);
            Cd = (ReD - 200000) / (1000000 - 200000) * (Cd2 - Cd1) + Cd1;
        } else if (ReD >= 10000000) { //REd is out of range of correlations
            Cd = 1 - 0.184 / pow(ReD, 0.2);
        } else { //REd is out of range of correlations
            Cd = 1 - 6.92 / sqrt(ReD);
        }
		return Cd;
}
