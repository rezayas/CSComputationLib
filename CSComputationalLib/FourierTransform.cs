using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using RandomNumberGeneratorLib;

namespace ComputationLib
{
    public static class FourierTransform
    {
        public static void DoFourierTransform(double[] data, ref double[] periods, ref double[] weights, ref double[] recon_timeseries)
        {
            periods = new double[(int)Math.Floor(data.Length / 2.0f) + 1];
            weights = new double[(int)Math.Floor(data.Length / 2.0f) + 1];
            recon_timeseries = new double[data.Length];

            alglib.complex[] cvalued_fft;
            alglib.complex[] orig_cvalued_fft;
            alglib.complex[] thresholded_cvalued_fft = new alglib.complex[(int)Math.Floor(data.Length / 2.0f) + 1];

            alglib.fftr1d(data, out orig_cvalued_fft);
            double[] amplitude = new double[(int)Math.Floor(data.Length / 2.0f) + 1];
            ComplexToAmplitude(orig_cvalued_fft, ref amplitude);

            //System.Random rnd = new System.Random();
            RandomNumberGeneratorLib.ThreadSpecificRNG myRND = new ThreadSpecificRNG();
            myRND.Reset(0);

            double[][] array2Da = new double[(int)Math.Floor(data.Length / 2.0f) + 1][];
            for (int j = 0; j < array2Da.GetLength(0); j++)
            {
                array2Da[j] = new double[1000];
            }

            for (int i = 0; i < 1000; i++)
            {
                double[] RandomPermuted_Input = data.OrderBy(x => myRND.RNDU01()).ToArray(); //rnd.Next()
                double[] tmp_amplitude = new double[(int)Math.Floor(data.Length / 2.0f) + 1];

                alglib.fftr1d(RandomPermuted_Input, out cvalued_fft);
                ComplexToAmplitude(cvalued_fft, ref tmp_amplitude);

                for (int j = 0; j < tmp_amplitude.Length; j++)
                    array2Da[j][i] = tmp_amplitude[j];
            }


            for (int j = 0; j < array2Da.GetLength(0); j++)
            {
                Array.Sort(array2Da[j]); // ascending
                if (array2Da[j][(int)Math.Floor(0.99f * array2Da[j].Length)] < amplitude[j])
                {
                    weights[j] = amplitude[j];
                    thresholded_cvalued_fft[j] = orig_cvalued_fft[j];
                }
                else
                {
                    weights[j] = 0;
                    thresholded_cvalued_fft[j] = 0;
                }
                if (j != 0)
                {
                    periods[j] = 2.0f * ((double)array2Da.GetLength(0)-1) / j;
                }
                else
                {
                    thresholded_cvalued_fft[j] = orig_cvalued_fft[j];
                }

            }

            alglib.fftr1dinv(thresholded_cvalued_fft, recon_timeseries.Length, out recon_timeseries);
        }

        static void ComplexToAmplitude(alglib.complex[] data, ref double[] amplitude)
        {
            for (int i = 0; i < (int)Math.Floor(data.Length / 2.0f) + 1; i++)
            {
                amplitude[i] = 2 * Math.Sqrt(Math.Pow(data[i].x, 2) + Math.Pow(data[i].y, 2));
            }
        }
    }
}
