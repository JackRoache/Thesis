__kernel
void besslj0(__global float *in, __global float *out, const int num)
{
    const int id = get_global_id(0);

    if (id < num){
        float x = in[id];
        float ax, z;
        float xx, y, ans, ans1, ans2, t;

        if ((ax = fabs(x)) < 8.0) {
            y = x*x;
            ans1 = 57568490574.0 + y*(-13362590354.0 + y*(651619640.7
                                                          + y*(-11214424.18 + y*(77392.33017 + y*(-184.9052456)))));
            ans2 = 57568490411.0 + y*(1029532985.0 + y*(9494680.718
                                                        + y*(59272.64853 + y*(267.8532712 + y*1.0))));
            ans = ans1 / ans2;

        } else {
            z = 8.0 / ax;
            y = z*z;
            xx = ax - 0.785398164;
            ans1 = 1.0 + y*(-0.1098628627e-2 + y*(0.2734510407e-4 + y*(-0.2073370639e-5 + y*0.2093887211e-6)));
            ans2 = -0.1562499995e-1 + y*(0.1430488765e-3 + y*(-0.6911147651e-5 + y*(0.7621095161e-6 - y*0.934935152e-7)));
            t = 0.636619772 / ax;
            ans = sqrt(t) *(cos(xx)*ans1 - z*sin(xx)*ans2);
        }
        out[id] = ans;
    }
}

__kernel
void besslj1(__global *in, __global float *out, const int num)
{
    const int id = get_global_id(0);

    if (id < num) {
        float x = in[id];
        float ax, z, t;
        float xx, y, ans, ans1, ans2;

        if ((ax = fabs(x)) < 8.0) {
            y = x*x;
            ans1 = x*(72362614232.0 + y*(-7895059235.0 + y*(242396853.1
                                                            + y*(-2972611.439 + y*(15704.48260 + y*(-30.16036606))))));
            ans2 = 144725228442.0 + y*(2300535178.0 + y*(18583304.74
                                                         + y*(99447.43394 + y*(376.9991397 + y*1.0))));
            ans = ans1 / ans2;
        }
        else {
            z = 8.0 / ax;
            y = z*z;
            xx = ax - 2.356194491;
            ans1 = 1.0 + y*(0.183105e-2 + y*(-0.3516396496e-4
                                             + y*(0.2457520174e-5 + y*(-0.240337019e-6))));
            ans2 = 0.04687499995 + y*(-0.2002690873e-3
                                      + y*(0.8449199096e-5 + y*(-0.88228987e-6
                                                                + y*0.105787412e-6)));
            t = 0.636619772 / ax;
            ans = sqrt(t)*(cos(xx)*ans1 - z*sin(xx)*ans2);
            if (x < 0.0) ans = -ans;
        }
        out[id] = ans;
    }
}

__kernel
void bessly0(__global float *in, __global float *out, const int num)
{
    const int id = get_global_id(0);

    if (id < num){
        float x = in[id];
        double z;
        double xx, y, ans, ans1, ans2, t;

        if (x < 8.0) {
                y = x*x;
                ans1 = -2957821389.0 + y*(7062834065.0 + y*(-512359803.6
                        + y*(10879881.29 + y*(-86327.92757 + y*228.4622733))));
                ans2 = 40076544269.0 + y*(745249964.8 + y*(7189466.438
                        + y*(47447.26470 + y*(226.1030244 + y*1.0))));
                ans = (ans1 / ans2) + 0.636619772*bessj0(x)*log(x);
        }
        else {
                z = 8.0 / x;
                y = z*z;
                xx = x - 0.785398164;
                ans1 = 1.0 + y*(-0.1098628627e-2 + y*(0.2734510407e-4
                        + y*(-0.2073370639e-5 + y*0.2093887211e-6)));
                ans2 = -0.1562499995e-1 + y*(0.1430488765e-3
                        + y*(-0.6911147651e-5 + y*(0.7621095161e-6
                                + y*(-0.934945152e-7))));
                t = 0.636619772 / x;
                ans = sqrt(t)*(sin(xx)*ans1 + z*cos(xx)*ans2);
        }
        out[id] = ans;
    }
}
