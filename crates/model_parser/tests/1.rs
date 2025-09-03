use model_parser::Builder;

#[derive(Builder, Debug)]
struct A {
    a: String,
    b: i32,
}

#[test]
fn test_a() {
    let a: A = A::builder().a("Hello".to_string()).b(42).build();
    println!("{:?}", a);
}
