//#  simulate samples from student t distribution 
//# 
// http://www.openbugs.net/Examples/t-df.html

transformed data {
  int d;
  d = 4;
}
parameters {
  array[1000] real y;
}
model {
  y ~ student_t(d, 0, 1);
}
