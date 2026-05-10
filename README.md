# Inferencia verificable en modelos de Aprendizaje Automático

**Verifiable Inference in Machine Learning Models**

> Trabajo académico para la asignatura de **Criptografía** — Escuela Superior de Ingeniería Informática, Universidad de Castilla-La Mancha (Albacete)
>
> Academic project for the **Cryptography** course — Escuela Superior de Ingeniería Informática, Universidad de Castilla-La Mancha (Albacete, Spain)

**Autor / Author:** Álvaro Martínez Alfaro (`alvaro.martinez39@alu.uclm.es`)

---

## Español

### Descripción

Este trabajo estudia la intersección entre el **aprendizaje automático** y las **pruebas de conocimiento cero** (ZKPs), campo conocido como *Zero-Knowledge Machine Learning* (ZKML). El punto de partida es el modelo de negocio de inferencia como servicio: un proveedor despliega un modelo entrenado y responde a consultas sin revelar sus pesos. Esto introduce un problema de confianza — el usuario no puede saber si el proveedor ejecutó realmente el modelo que afirma ejecutar, ni si lo hizo correctamente — mientras que el proveedor tiene incentivos legítimos para mantener sus pesos en secreto. Las pruebas de conocimiento cero permiten demostrar que un cómputo se ha realizado correctamente sin revelar información secreta sobre él, ofreciendo así una solución teórica a esta tensión.

El trabajo cubre:

- **Fundamentos de redes neuronales artificiales**: inferencia, entrenamiento y funciones de activación.
- **Pruebas de conocimiento cero**: pruebas interactivas, sistemas de prueba, SNARKs y zk-SNARKs.
- **ZKML — desafíos técnicos**:
  - Aritmetización de la inferencia mediante sistemas **R1CS** y programas aritméticos cuadráticos (**QAP**).
  - **Cuantización** de pesos y activaciones a aritmética de cuerpo finito.
  - Manejo de funciones de activación no lineales (e.g., **ReLU**).
- **Revisión de frameworks modernos**: EZKL, Circom, ZEN y zkLLM.
- **Prueba de concepto en MATLAB**: implementación didáctica del flujo completo (cuantización → R1CS → QAP → verificación polinómica).

### Estructura del repositorio

| Fichero / Directorio | Descripción |
|---|---|
| `main.tex` | Documento LaTeX principal de la memoria |
| `main.pdf` | Versión compilada de la memoria |
| `refs.bib` | Referencias bibliográficas en BibTeX |
| `zkml_poc.mlx` | Implementación en MATLAB (*Live Script*) |
| `zkml_poc.pdf` | Exportación PDF de la implementación MATLAB |
| `figs/` | Figuras e imágenes utilizadas en la memoria |

### Prueba de concepto en MATLAB

La implementación (`zkml_poc.mlx`) reproduce de forma simplificada el proceso completo de inferencia verificable sobre una red neuronal pequeña (arquitectura 2×2×1, activación ReLU, cuerpo finito F₆₅₅₃₇):

1. **Cuantización**: conversión de pesos y activaciones a aritmética modular.
2. **Aritmetización R1CS**: traducción de cada operación de la inferencia a restricciones R1CS.
3. **Transformación a QAP**: obtención de los polinomios A, B, C mediante interpolación de Lagrange.
4. **Verificación**: generación del polinomio cociente H y verificación de la identidad polinómica por evaluación en un punto aleatorio.

> Nota: la implementación constituye un SNARK sobre la inferencia, pero no un zk-SNARK completo, ya que los componentes que aportan la propiedad de conocimiento cero (compromisos criptográficos, grupos de orden primo grande, emparejamientos bilineales) quedan fuera del alcance didáctico.

### Licencia

El código fuente está licenciado bajo MIT (ver `LICENSE`).
La memoria está licenciada bajo 
[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).

---

## English

### Description

This project studies the intersection of **machine learning** and **zero-knowledge proofs** (ZKPs), a field known as *Zero-Knowledge Machine Learning* (ZKML). The starting point is the inference-as-a-service business model: a provider deploys a trained model and responds to queries without revealing its weights. This creates a trust problem — the user has no way to know whether the provider actually ran the model it claims to run, or whether it did so correctly — while the provider has legitimate incentives to keep its weights secret. Zero-knowledge proofs allow one party to convince another that a computation was performed correctly without revealing any secret information about it, offering a theoretical solution to this tension.

The report covers:

- **Artificial neural network fundamentals**: inference, training, and activation functions.
- **Zero-knowledge proofs**: interactive proofs, proof systems, SNARKs, and zk-SNARKs.
- **ZKML — technical challenges**:
  - Arithmetization of inference via **R1CS** constraint systems and **Quadratic Arithmetic Programs (QAP)**.
  - **Quantization** of weights and activations to finite-field arithmetic.
  - Handling non-linear activation functions (e.g., **ReLU**) via bit decomposition and lookup tables.
- **Review of modern frameworks**: EZKL, Circom, ZEN, and zkLLM.
- **MATLAB proof of concept**: a didactic implementation of the full pipeline (quantization → R1CS → QAP → polynomial verification).

### Repository structure

| File / Directory | Description |
|---|---|
| `main.tex` | Main LaTeX source of the report |
| `main.pdf` | Compiled PDF of the report |
| `refs.bib` | BibTeX bibliography |
| `zkml_poc.mlx` | MATLAB implementation (*Live Script*) |
| `zkml_poc.pdf` | PDF export of the MATLAB implementation |
| `figs/` | Figures and images used in the report |

### MATLAB proof of concept

The implementation (`zkml_poc.mlx`) reproduces, in a simplified form, the full verifiable inference pipeline over a small neural network (2×2×1 architecture, ReLU activation, finite field F₆₅₅₃₇):

1. **Quantization**: convert weights and activations to modular arithmetic.
2. **R1CS arithmetization**: translate each inference operation into R1CS constraints.
3. **QAP transformation**: obtain polynomials A, B, C via Lagrange interpolation.
4. **Verification**: compute the quotient polynomial H and verify the polynomial identity by evaluation at a random point.

> Note: the implementation constitutes a SNARK over the inference, but not a full zk-SNARK, since the components that provide the zero-knowledge property (cryptographic commitments, large-prime-order groups, bilinear pairings) are out of scope for this didactic exercise.

### License

The source code is distributed under the MIT License (see `LICENSE`).
The accompanying report is released under 
[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).
