use rand::Rng;
const PI: f64 = std::f64::consts::PI;

type KuramotoArgs = (Vec<f64>, f64);

fn cauchy(gamma: f64) -> f64 {
    let mut rng = rand::thread_rng();
    let x: f64 = rng.gen();
    let y: f64 = gamma*(x-0.5).tan();
    return y;
}

fn theta_init() -> f64 {
    let mut rng = rand::thread_rng();
    let x: f64 = rng.gen();
    let y: f64 = 2.0*PI*x;
    return y;
}

//vector field of kuramoto model
fn kuramoto(xs: &Vec<f64>, args: &KuramotoArgs) -> Vec<f64>{
    let omegas = &args.0;
    let K = &args.1;
    let mut vs: Vec<f64> = omegas.clone();
    let mut rsin: f64 = 0.0;
    let mut rcos: f64 = 0.0;
    for x in xs.iter(){
        rsin+=x.sin();
        rcos+=x.cos();
    }
    rsin/=xs.len() as f64;
    rcos/=xs.len() as f64;
    for (i, x) in xs.iter().enumerate(){
        vs[i] += K*rsin*x.cos() - K*rcos*x.sin();
    }
    return vs;
}

fn add_vec(a: &Vec<f64>, b: &Vec<f64>) -> Vec<f64>{
    let mut vs = Vec::with_capacity(a.len());
    for i in 0..a.len(){
        vs.push(a[i] + b[i]);
    }
    return vs;
}

fn multiple_vec(a: &Vec<f64>, c: f64) -> Vec<f64>{
    let mut vs = Vec::with_capacity(a.len());
    for i in 0..a.len(){
        vs.push(a[i] * c);
    }
    return vs;
}

fn order_parameter(thetas: &Vec<f64>) -> f64{
    let mut rsin: f64 = 0.0;
    let mut rcos: f64 = 0.0;
    for t in thetas.iter(){
        rsin+=t.sin();
        rcos+=t.cos();
    }
    rsin/=thetas.len() as f64;
    rcos/=thetas.len() as f64;
    return (rsin*rsin+rcos*rcos).sqrt();
}

//Argsはvector_fieldの引数として取る
fn runge_kutta<Args>(vector_fields: fn(&Vec<f64>, &Args) -> Vec<f64>, xs: Vec<f64>, args: &Args, dt: f64) -> Vec<f64> {
    let k1 = vector_fields(&xs, &args);
    let k2 = vector_fields(&add_vec(&xs, &multiple_vec(&k1, dt*0.5)), &args);
    let k3 = vector_fields(&add_vec(&xs, &multiple_vec(&k2, dt*0.5)),&args);
    let k4 = vector_fields(&add_vec(&xs, &multiple_vec(&k3, dt)), &args);
    // calc slope from k1, k2, k3 and k4
    let slope = multiple_vec(&add_vec(&add_vec(&k1, &k4), &multiple_vec(&add_vec(&k2, &k3), 2.0)), dt/6.0);
    // calc next step
    let next_xs = add_vec(&xs, &slope);
    return next_xs;
}

fn main() {
    //sizeは将来はコマンドライン引数にする
    const SIZE: usize = 10000;
    let gamma: f64 = 1.0;
    let mut omegas: Vec<f64> = Vec::with_capacity(SIZE);
    let mut thetas: Vec<f64> = Vec::with_capacity(SIZE);
    for _i in 0..SIZE{
        omegas.push(cauchy(gamma));
        thetas.push(theta_init());
    }
    let K: f64 = 0.6;
    let k_a: KuramotoArgs = (omegas, K);
    let dt: f64 = 0.05;
    let loops: i64 = 2000;
    let mut rs: Vec<f64> = Vec::with_capacity(loops as usize);
    
    for _i in 0..loops {
        thetas = runge_kutta::<KuramotoArgs>(kuramoto, thetas, &k_a, dt);
        rs.push(order_parameter(&thetas));
    }
    
    println!("Hello, world");
    
    println!("{:?}", rs);
}