t1 = cos(xi);
t2 = cos(eta);
t3 = cos(K1bry);
t4 = sin(xi);
t5 = sin(K1bry);
t6 = sin(K1brx);
t7 = sin(eta);
t8 = cos(K1brx);
t9 = t2 * t4;
t10 = -t9 * t6 - t7 * t8;
t11 = t1 * t2;
t12 = sin(zeta);
t13 = cos(zeta);
t14 = t1 * t12;
t15 = t4 * t13;
t16 = t14 * t7 + t15;
t4 = t4 * t12;
t1 = t1 * t13;
t17 = -t4 * t7 + t1;
t18 = t2 * t8;
t19 = t18 * t12 + t17 * t6;
t1 = -t1 * t7 + t4;
t4 = t15 * t7 + t14;
t14 = -t18 * t13 + t4 * t6;
cg0 = [t10 * t5 + t11 * t3 t16 * t3 + t19 * t5 t1 * t3 + t14 * t5; t6 * t7 - t9 * t8 -t2 * t12 * t6 + t17 * t8 t2 * t13 * t6 + t4 * t8; -t10 * t3 + t11 * t5 t16 * t5 - t19 * t3 t1 * t5 - t14 * t3;];