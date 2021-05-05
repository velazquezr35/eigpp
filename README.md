# ayudantia_eigpp_2021


### Apéndice original del código

La clase principal y sus dos subclases asociadas se definen por:

- sim:
-     stru
-     aero

Reservando sendos casos para análisis estructurales y con cargas aerodinámicas, respectivamente.

La información se guarda bajo las siguientes variables, pertenecientes a la subclase estructural:

- sim.stru.mass: Reales, ngl - matriz de masa diagonal del modelo.
- sim.stru.om: Reales, nm - vector de frecuencias naturales ordenadas de menor a mayor.
- sim.stru.phi: Reales, ngl x nm - matriz modal (cada columna es un modo).
- sim.stru.q: Reales, ngl x nt - matriz con coordenadas modales a lo largo del tiempo, cada columna tiene el vector q de un instante.
- sim.stru.t: Reales, nt - vector con instantes donde se conoce la solución (malla o grilla temporal).
- sim.stru.nt: Entero, 1 - cantidad de instantes de tiempo.
- sim.stru.u raw: Reales, ngl x nt - matriz con desplazamientos a lo largo del tiempo, el vector u para un instante se guarda en una columna - como se obtienen de Simpact/curvas (si existen GL rotacionales, se expresan con ángulos de Euler y no se puede usar para hacer descomposición modal).
- sim.stru.u avr: Reales, ngl x nt - DyG con GR convertidos a expresión vectorial.
- sim.name: Nombre de la simulación (se usa luego para guardar archivos y demás).
- sim.stru.nodes: Etiquetas de nodos de interés - lista de enteros.
- sim.stru.nnode: Cantidad de nodos considerados.
- sim.stru.ndof: Cantidad de GL considerados (incluyendo los que son restringidos por condiciones de borde).
- sim.stru.rdOpt: Bandera para leer datos.
- sim.stru.rdof: Bandera que indique si hay GL rotacionales.
- sim.str.eigOpt: Bandera para hacer (o no) descomposición modal.
