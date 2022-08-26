using Pumas # branch mt/corr_cholesky2
using DataFrames

poppk2cpt = @model begin
  @param begin
    tvcl ~ LogNormal(log(10), 0.25) # CL
    tvq ~ LogNormal(log(15), 0.5)   # Q
    tvvc ~ LogNormal(log(35), 0.25) # V1
    tvvp ~ LogNormal(log(105), 0.5) # V2
    tvka ~ LogNormal(log(2.5), 1)   # ka
    σ ~ truncated(Cauchy(0, 5), 0, Inf) # sigma
    C ~ LKJCholesky(5, 1.0)
    ω ∈ VectorDomain(5; lower=0.01, upper=2.0)
  end

  @random begin
    ηstd ~ MvNormal(I(5))
  end

  @pre begin
    η = ω .* (getchol(C).L * ηstd)
    CL = tvcl * exp(η[1])
    Q = tvq * exp(η[2])
    Vc = tvvc * exp(η[3])
    Vp = tvvp * exp(η[4])
    Ka = tvka * exp(η[5])
  end

  @dynamics Depots1Central1Periph1

  @derived begin
    # Torsten uses log(dv) ~ Normal(log(cp), sigma)
    cp := @. Central / Vc
    dv ~ @. LogNormal(log(cp), σ)
  end
end

# Torsten Data
addl = [14, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 14, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 14, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  14, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  14, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  14, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  14, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  14, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  14, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  14, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
amt = [80000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 80000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 80000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 80000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 80000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 80000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 80000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 80000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 80000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 80000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0]
cmt = [1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
  2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
  2, 2, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
  2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
  1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
  2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 2,
  2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
  2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 2, 2, 2,
  2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
  2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 2, 2, 2, 2, 2,
  2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
  2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 2, 2, 2, 2, 2, 2, 2,
  2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
  2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2,
  2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
  2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
  2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
  2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
  2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
  2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]
cObs = [
  missing,
  397.975649907265, 662.021672889522, 1034.24388828986, 1450.22976548349,
  1999.41045079086, 1866.34871848592, 1817.8552679785, 1387.28649787039,
  1177.05185401853, 786.069845789358, 468.879464944515, 438.973120196479,
  266.626209140417, 632.998498582195, 842.544242531116, 1131.60117022217,
  1756.13450408554, 2425.3700813798, 2302.24085807503, 2364.24277477216,
  2110.02937116503, 1344.94459907031, 970.361876344979, 755.4020091503,
  518.250815167332, 444.03875806556, 630.430411973066, 571.198353856031,
  645.887738000823, 773.34397946794, 740.86621556372, 653.351484128769,
  733.419269348437, 822.018722686982, 786.120707361446, 829.310743285099,
  818.041014594213, 747.237833718697, 1345.52700130719, 1796.87452544178,
  1761.53332476314, 2701.86366864204, 3276.63831445814, 2721.69796730839,
  2428.73711355341, 2619.82777613324, 1666.14062302942, 1438.98037540228,
  813.467412780766, 914.22254128572, 673.756833007774, 529.725577945499,
  483.562385175449,
  missing,
  436.197093426514, 611.279654955084, 920.908255945113,
  1627.08828521263, 1786.2553083795, 1738.81714746253, 1455.33769391704,
  1692.71689547075, 1255.85237954741, 885.789601059148, 448.360127274659,
  363.377355043768, 319.263045347803, 585.88040406814, 996.195737398005,
  1329.59040223411, 1721.87295474567, 2105.68744489902, 2292.05632341759,
  2061.62129009226, 1881.65597612589, 1426.61745339035, 1023.60796451466,
  770.788025273378, 572.199647515165, 473.66026623036, 562.672189089989,
  651.689694716912, 666.490458293232, 649.051236845097, 730.792276390941,
  881.995031900801, 657.065288611115, 744.828213075366, 723.20788579624,
  749.414465493754, 675.733863036767, 769.703519350558, 1139.8899387404,
  1426.84498363203, 1859.11207314864, 2665.30097882561, 2701.23294190221,
  2706.42537556177, 2732.153131495, 2811.14599022345, 1714.80037516719,
  1204.14533481773, 1167.53823934818, 948.610044150629, 797.427352166287,
  655.854972196013, 463.653624343523,
  missing,
  367.817745204874, 774.22374633784,
  1067.22283253391, 1676.05538128583, 1697.67178988455, 2121.63227045589,
  2072.36383197827, 1654.62030314725, 1155.24340944936, 693.744644005515,
  429.988549712065, 291.130916625659, 260.434605759822, 547.748071792611,
  1030.38076228306, 1215.78640585911, 2076.42386410568, 2534.65698204429,
  2255.12436598015, 2116.96950824167, 1844.31641358679, 1456.89893920159,
  981.555036143453, 695.661606524004, 548.31127110545, 448.521915325858,
  569.44570301264, 649.172883065837, 631.666949583856, 640.464895472125,
  665.871284410796, 796.394698865376, 733.706946267985, 680.179857126242,
  870.751491966711, 861.502325886929, 758.386243709701, 777.513507169934,
  1153.67155184892, 1242.36603116282, 1822.5758913193, 2292.01656553055,
  2843.01051401753, 2518.55458802068, 2971.283561892, 2685.73717287768,
  1843.78933394958, 1318.37278980453, 1338.57035730673, 883.712083980592,
  800.874783338039, 642.444142666784, 479.616109763875,
  missing,
  354.953883930517,
  763.514916353009, 829.171504748064, 1461.35190200939, 1826.34415694983,
  2008.26137002812, 1845.84318429145, 1904.08270350537, 1204.24628198449,
  846.943652649601, 445.202319720968, 387.084799629131, 310.067190039326,
  581.539328627347, 1138.18539341057, 1432.58984467926, 1751.22596150434,
  2148.3510876309, 2195.73418672597, 1796.06233461043, 2165.10917101154,
  1289.88101744447, 1170.18206921783, 676.539247805496, 548.854737715062,
  429.558563143288, 574.370704571517, 740.560521599011, 470.993415248098,
  708.583400726987, 823.172319541011, 770.112562351309, 733.96490993191,
  687.679795553633, 726.137880443138, 800.253282924768, 734.633120687101,
  812.584319181785, 1168.49296872008, 1320.16365121206, 1756.50721828332,
  2025.34509141398, 2527.71259431052, 2926.23466943579, 2733.21185537613,
  2361.85449720578, 1981.56687238079, 1382.15882782796, 1043.09628490683,
  915.509131365557, 897.41695460627, 644.519341459823, 448.62881278082,
  missing,
  353.009509314532, 789.298193947955, 841.325050084809, 1558.90124992445,
  2196.3612760194, 1910.13893840071, 2383.88822698888, 1611.87283945709,
  1127.51205814178, 906.159224705914, 428.912506016965, 392.293631064577,
  283.074811595929, 721.294501135592, 1170.526253347, 977.329390137667,
  1697.23594271301, 1999.12997790124, 1984.09917693634, 1902.34821449105,
  1922.14294569924, 1439.68971955885, 1073.02995920878, 730.469424479868,
  573.336208734975, 413.713167594117, 454.31954652466, 666.039057042795,
  759.053424259795, 658.014657573545, 832.399304836682, 793.635989116005,
  788.00395813532, 758.503979651513, 819.499225815918, 725.588178447319,
  607.562307221738, 659.817967136868, 1123.95163614892, 1579.00209018514,
  1510.7440418548, 2090.26444121586, 2625.04813504285, 2908.88006551352,
  2601.5241573885, 2580.74672665599, 1772.30045871553, 1573.91193237866,
  1076.84470265722, 852.250795614213, 755.213353302628, 666.816672598448,
  494.095545149186,
  missing,
  353.217833756205, 738.741480519181, 892.979770247455,
  1646.50898566913, 1803.65145489691, 1844.62585852542, 2052.19137089459,
  1486.75144532293, 1179.89799503594, 758.038613010266, 442.194322320324,
  370.269457105073, 258.614266165473, 660.084598190974, 951.15594669722,
  1081.82356095063, 1897.65092925175, 2126.0198062466, 1707.24612330576,
  2273.72699999375, 1897.99180625451, 1355.1239596411, 1133.38375334675,
  703.900735726078, 569.610977095857, 524.313640177172, 593.159471560388,
  620.875892680826, 673.355964796036, 664.863955932219, 786.657596557657,
  799.198531637453, 717.573071519453, 751.427936834861, 782.549544689467,
  803.765333726231, 703.548945115896, 828.98557490575, 932.295007241806,
  1553.35529610592, 1589.56286584338, 2152.06794498527, 2361.64263469279,
  2468.96238077794, 2717.49547599601, 2278.06184466794, 1688.71969106609,
  1339.95441355306, 849.426694403967, 1022.4022680086, 740.225261597107,
  659.909668740848, 471.481608111606,
  missing,
  398.548240800401, 727.267246547714,
  964.125363045441, 1491.58206980841, 1752.69301224089, 1996.87360732426,
  1765.06365935688, 1528.89121192213, 1064.61097922577, 781.964631776258,
  490.343906402351, 353.14891560134, 217.748305527967, 586.954027953212,
  938.264614484194, 1301.56581745975, 1656.59017083135, 1784.55091639707,
  1873.28207598533, 2196.70092996848, 1844.74262704192, 1410.81295081477,
  984.511392134695, 691.658035282262, 615.606564715515, 474.915680896039,
  525.101177925788, 549.264755302626, 673.82967191078, 689.454392161415,
  733.192493666863, 693.476765539752, 728.639631833443, 669.635244930601,
  788.592285370643, 659.783770272102, 732.25392435712, 779.457041433298,
  982.509269952802, 1253.1073008132, 1914.83980785555, 2241.07947102336,
  2272.68081142453, 2918.81292573619, 3174.91284741816, 2230.28033644049,
  1672.95004886109, 1612.63468017055, 1192.66521461374, 942.061951945441,
  807.445040328342, 616.372550526072, 514.975941088068,
  missing,
  345.831597652523,
  706.798564928919, 998.460548053518, 1613.51076611719, 2452.22487509151,
  2172.69069981241, 1700.88326723288, 1589.39668967283, 1133.53715224905,
  729.102983379094, 532.478370179644, 374.19541855252, 278.759453606893,
  720.668659266523, 980.738111970951, 1113.88197839311, 1471.84530425725,
  2310.14876947441, 2300.64847767868, 1765.56776792361, 1767.52763359072,
  1325.08968391444, 948.914315361711, 878.959502496669, 476.560218512476,
  382.034581491127, 629.865230521741, 685.517839492663, 607.634552490588,
  640.081110357082, 691.603843420679, 698.089005683888, 647.035076360675,
  707.026578055257, 749.596141930849, 644.135131004664, 761.360130349336,
  766.469990509427, 1002.94212352628, 1338.44300357825, 1673.04256212687,
  2824.23654374767, 2442.28846014567, 3037.95268117384, 2462.39922643292,
  2562.61372616933, 2253.66360752813, 1561.99757035381, 872.179550789084,
  926.244154115525, 688.796845995604, 559.114772264193, 531.612085580992,
  missing,
  389.702733762232, 673.059666771206, 808.453597577418, 1571.80544976789,
  1686.04216700093, 1669.81062580809, 2080.06077511601, 1941.43545624699,
  1013.86969305464, 944.120766167887, 585.366077149881, 346.857395948262,
  259.930103228702, 564.976770859736, 980.323010910724, 1144.16901704834,
  1727.38496006909, 2193.79320885309, 2225.65754022282, 1743.27137490536,
  1964.12493380825, 1289.34110641288, 936.712663151495, 626.747213867843,
  723.566892124005, 396.160622594282, 508.202083498349, 666.420141035693,
  717.576053730096, 739.105368107922, 698.320931972822, 797.189502479437,
  622.313920726652, 742.248821015959, 688.404476738483, 776.316632419954,
  741.409648737448, 682.605264190909, 1290.85801022963, 1774.0905791404,
  1760.63747852042, 2124.71973884039, 2615.32736325754, 2750.55716568369,
  2947.15142449091, 3080.74272117931, 1699.95918088019, 1620.89875779404,
  978.015370598894, 936.8038279724, 747.011522101322, 560.50782459446,
  449.340015090543,
  missing,
  341.725376545919, 764.096248760819, 1069.63251859874,
  1781.72166816296, 1778.9491335381, 1714.87085587664, 1883.9886038498,
  1603.22016191613, 1158.94386028094, 861.401653595964, 491.449734996973,
  341.87762841343, 280.556076982054, 629.089217179725, 951.083010898348,
  1327.77337309185, 1943.63489952064, 2324.63200035449, 2441.31856277739,
  2418.36876841171, 1864.22181704412, 1329.89709285803, 1005.66707180987,
  721.689631542447, 576.98909530742, 427.627774781389, 628.361002142369,
  650.073117805091, 664.260294701744, 653.326862506559, 728.24773306327,
  695.855593291352, 727.582614532649, 749.917551293357, 752.44981928654,
  677.627083606345, 793.54857332707, 774.784886290835, 868.683055860495,
  1503.33067922814, 1714.82663677575, 2034.91802344815, 2444.04293660476,
  3150.77258725255, 2766.87887817365, 2129.02524750355, 1771.01585543452,
  1277.03135200084, 951.330480187572, 1019.59765492059, 725.076650100637,
  641.002132021179, 406.047208562917]
evid = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
ii = [12, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 12, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 12, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  12, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  12, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  12, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  12, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  12, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  12, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  12, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
nSubjects = 10
rate = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
ss = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
time = [0, 0.083, 0.167, 0.25, 0.5, 0.75, 1, 1.5, 2, 3, 4, 6, 8, 12, 12.083,
  12.167, 12.25, 12.5, 12.75, 13, 13.5, 14, 15, 16, 18, 20, 24, 36, 48, 60, 72, 84,
  96, 108, 120, 132, 144, 156, 168, 168.083, 168.167, 168.25, 168.5, 168.75, 169,
  169.5, 170, 171, 172, 174, 176, 180, 186, 192, 0, 0.083, 0.167, 0.25, 0.5, 0.75, 1,
  1.5, 2, 3, 4, 6, 8, 12, 12.083, 12.167, 12.25, 12.5, 12.75, 13, 13.5, 14, 15, 16,
  18, 20, 24, 36, 48, 60, 72, 84, 96, 108, 120, 132, 144, 156, 168, 168.083, 168.167,
  168.25, 168.5, 168.75, 169, 169.5, 170, 171, 172, 174, 176, 180, 186, 192, 0,
  0.083, 0.167, 0.25, 0.5, 0.75, 1, 1.5, 2, 3, 4, 6, 8, 12, 12.083, 12.167, 12.25,
  12.5, 12.75, 13, 13.5, 14, 15, 16, 18, 20, 24, 36, 48, 60, 72, 84, 96, 108, 120,
  132, 144, 156, 168, 168.083, 168.167, 168.25, 168.5, 168.75, 169, 169.5, 170, 171,
  172, 174, 176, 180, 186, 192, 0, 0.083, 0.167, 0.25, 0.5, 0.75, 1, 1.5, 2, 3, 4, 6,
  8, 12, 12.083, 12.167, 12.25, 12.5, 12.75, 13, 13.5, 14, 15, 16, 18, 20, 24, 36,
  48, 60, 72, 84, 96, 108, 120, 132, 144, 156, 168, 168.083, 168.167, 168.25, 168.5,
  168.75, 169, 169.5, 170, 171, 172, 174, 176, 180, 186, 192, 0, 0.083, 0.167, 0.25,
  0.5, 0.75, 1, 1.5, 2, 3, 4, 6, 8, 12, 12.083, 12.167, 12.25, 12.5, 12.75, 13, 13.5,
  14, 15, 16, 18, 20, 24, 36, 48, 60, 72, 84, 96, 108, 120, 132, 144, 156, 168,
  168.083, 168.167, 168.25, 168.5, 168.75, 169, 169.5, 170, 171, 172, 174, 176, 180,
  186, 192, 0, 0.083, 0.167, 0.25, 0.5, 0.75, 1, 1.5, 2, 3, 4, 6, 8, 12, 12.083,
  12.167, 12.25, 12.5, 12.75, 13, 13.5, 14, 15, 16, 18, 20, 24, 36, 48, 60, 72, 84,
  96, 108, 120, 132, 144, 156, 168, 168.083, 168.167, 168.25, 168.5, 168.75, 169,
  169.5, 170, 171, 172, 174, 176, 180, 186, 192, 0, 0.083, 0.167, 0.25, 0.5, 0.75, 1,
  1.5, 2, 3, 4, 6, 8, 12, 12.083, 12.167, 12.25, 12.5, 12.75, 13, 13.5, 14, 15, 16,
  18, 20, 24, 36, 48, 60, 72, 84, 96, 108, 120, 132, 144, 156, 168, 168.083, 168.167,
  168.25, 168.5, 168.75, 169, 169.5, 170, 171, 172, 174, 176, 180, 186, 192, 0,
  0.083, 0.167, 0.25, 0.5, 0.75, 1, 1.5, 2, 3, 4, 6, 8, 12, 12.083, 12.167, 12.25,
  12.5, 12.75, 13, 13.5, 14, 15, 16, 18, 20, 24, 36, 48, 60, 72, 84, 96, 108, 120,
  132, 144, 156, 168, 168.083, 168.167, 168.25, 168.5, 168.75, 169, 169.5, 170, 171,
  172, 174, 176, 180, 186, 192, 0, 0.083, 0.167, 0.25, 0.5, 0.75, 1, 1.5, 2, 3, 4, 6,
  8, 12, 12.083, 12.167, 12.25, 12.5, 12.75, 13, 13.5, 14, 15, 16, 18, 20, 24, 36,
  48, 60, 72, 84, 96, 108, 120, 132, 144, 156, 168, 168.083, 168.167, 168.25, 168.5,
  168.75, 169, 169.5, 170, 171, 172, 174, 176, 180, 186, 192, 0, 0.083, 0.167, 0.25,
  0.5, 0.75, 1, 1.5, 2, 3, 4, 6, 8, 12, 12.083, 12.167, 12.25, 12.5, 12.75, 13, 13.5,
  14, 15, 16, 18, 20, 24, 36, 48, 60, 72, 84, 96, 108, 120, 132, 144, 156, 168,
  168.083, 168.167, 168.25, 168.5, 168.75, 169, 169.5, 170, 171, 172, 174, 176, 180,
  186, 192]

# sanity check because I've added `missing` by hand in the cObs
evid1s = findall(i -> i == 1, evid)
@assert all((ismissing(cObs[i]) for i in evid1s)) # this should pass

df = DataFrame(; time, id=repeat(1:nSubjects; inner=54), ii, evid, dv=cObs, cmt, addl, amt)

# Pumas Population
pop = read_pumas(df)

# Inits from Torsten
# TODO: After fix bug order inits same as Torsten
init_params = (;
  tvcl=9.03101310740243,
  tvq=14.0196613638163,
  tvvc=84.7122459708543,
  tvvp=78.8103049860171,
  tvka=1.00344648906815,
  σ=1.23593237181194,
  C=float.(Matrix(I(5))),
  ω=[0.802001132862642, 1.09233685361687, 1.7213969589211, 1.08648055442609, 0.694992449134588]
)

poppk2cpt_fit = fit(poppk2cpt,
  pop,
  init_params,
  Pumas.BayesMCMC(nsamples=2_000, nadapts=1_000, target_accept=0.8, nchains=4))

Pumas.truncate(poppk2cpt_fit; burnin=1_000)
