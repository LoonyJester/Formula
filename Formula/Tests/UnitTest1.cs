using System;
using Formula;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace Tests
{
    [TestClass]
    public class UnitTest1
    {
        [TestMethod]
        public void TestMethod1()
        {
            var calculator = new CpCalculator();
            

            string ep = calculator.EP(CpCalculator.ProcedureType.Depletion, CpCalculator.GenderType.Female,  120, CpCalculator.WeightType.lbs,60,
                CpCalculator.HeightType.inc, 1000, CpCalculator.BloodVolumeType.Nadler, 40, 29, 30, 55, 1500, CpCalculator.ReplacementVolumeType.Calculate,
                21, 40, 1);


            Assert.IsNotNull(ep);
             Assert.AreNotEqual(int.Parse(ep), 0);
        }
    }
}
