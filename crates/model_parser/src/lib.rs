use proc_macro::TokenStream;
use quote::{format_ident, quote};

#[proc_macro_derive(Builder, attributes(builder))]
pub fn derive(input: TokenStream) -> TokenStream {
    // get name of the struct
    let input = syn::parse_macro_input!(input as syn::DeriveInput);
    let ident = &input.ident;

    let builder_ident = format_ident!("{}Builder", input.ident);
    
    // Extract fields from the struct
    let fields = match input.data {
        syn::Data::Struct(data_struct) => {
            match data_struct.fields {
                syn::Fields::Named(fields_named) => fields_named.named,
                _ => panic!("Builder derive only supports structs with named fields"),
            }
        }
        _ => panic!("Builder derive only supports structs"),
    };
    
    // Generate builder methods for each field
    let funcs = fields.iter().map(|field| {
        let field_name = field.ident.as_ref().expect("Expected named field");
        let field_type = &field.ty;

        quote! {
            pub fn #field_name(mut self, value: impl Into<#field_type>) -> Self {
                self.#field_name = Some(value.into());
                self
            }
        } 
    });
    
    // Generate builder struct fields
    let builder_fields = fields.iter().map(|field| {
        let field_name = field.ident.as_ref().expect("Expected named field");
        let field_type = &field.ty;
        
        quote! {
            #field_name: Option<#field_type>
        }
    });
    
    // Generate build method field assignments
    let build_fields = fields.iter().map(|field| {
        let field_name = field.ident.as_ref().expect("Expected named field");
        
        quote! {
            #field_name: self.#field_name.expect(&format!("Field '{}' is required", stringify!(#field_name)))
        }
    });

    quote! {
        impl #ident {
            pub fn builder() -> #builder_ident {
                #builder_ident::default()
            }
        }

        #[derive(Default)]
        struct #builder_ident {
            #(#builder_fields,)*
        }
        
        impl #builder_ident {
            pub fn build(self) -> #ident {
                #ident {
                    #(#build_fields,)*
                }
            }
            #(#funcs)*
        }
    }
    .into()
}
