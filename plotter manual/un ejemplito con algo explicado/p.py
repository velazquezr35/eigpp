
# import numpy as np
from sim_db import sim, rd_u, rd_rsn_De, svBin, modalDecomp, rdBin, nodeDof2idxEng
from plotter import fig_general, fig_FFT


bin_path = 'binaries/'
    
# create binaries
if False:
    objPath='corrida/'
    es=sim()
    es.stru.loadEigOpt=False
    es.name='example simulation'
    es.fName='lin_noAero'
    es.stru.nodes=list(range(200000,200036))
    es.stru.intLabOffset = 74
    es.stru.p11FN = 'snl.@1' #binary .p11 fname
    es.stru = rd_u(es.stru, **{'subDir_P11':objPath})
    es.stru.rsnDe= 'snl'
    es.stru.mnorm='norm' #'stiff'
    es.stru = rd_rsn_De(es.stru, **{'subDir_RSN':'modos/'})
    svBin(es, **{'subDir_BIN': bin_path})

if True:
    
    '''
    ALGUNOS COMENTARIOS:
        fig_general -----
            Necesita al menos tres argumentos:
            - un objeto de clase stru
            - una lista de listas con índices para buscar info
            - una lista análoga a la de índices, pero de atributos del objeto stru, dados como string,
                para plotear eso, dado como key-worded argument "data_type".
            Es decir, si quiero plotear la evolución temporal del modo 2 debo llamar así:
                fig_general(es.stru, [[2]], **{'data_type':[['q']]}).
            
            ATENCIÓN:
                notar que Python indexa desde 0 y, sin embargo, en la llamada a la función usé 2 para el segundo modo;
                esto es porque fig_general le resta uno a los índices para que la llamada tenga más sentido ingenieril.
            
            Que el segundo y el tercer argumentos sean "listas de listas" permite ordenar varias curvas en una misma
            figura de diferentes maneras. Por ejemplo, si quiero graficar el segundo modo en un par de ejes y
            el tercer modo en otro sistema de ejes, apilando verticalmente los dos gráficos, debo llamar como:
                fig_general(es.stru, [[2],[3]], **{'data_type':[['q'],['q']]});
            pero si quisiera graficarlos ambos en el mismo sistema de ejes debería usar:
                fig_general(es.stru, [[2,3]], **{'data_type':[['q','q']]});
            y, finalmente, si quiesiera agregar a la última versión dos ejes con los modos 10 y 12 debería hacer:
                fig_general(es.stru, [[2,3],[10],[12]], **{'data_type':[['q','q'],['q'],['q']]}).
            
            Los modos y todo lo afín a ellos se indexan de manera simple, de 1 al máximo, sin crear confusiones.
            Sin embargo, los GLG (grados de libertad geométricos) y todo lo afín, se indexan más naturalmente 
            considerando el número de nudo y luego el GL con una numeración local; por ejemplo, el desplazamiento
            vertical de la puntera es el GL 3 del nodo 200035. Este resulta ser el GL 213 de todo el sistema
            (esta es la numeración global).
            Existe una función que calcula (automáticamente) los índices de cada GLG a partir su nodo y su índice
            local: nodeDof2idxEng. Como argumento se necesitan un objeto de clase stru y un diccionario de la
            forma {nodo_A:[índices locales], nodo_B:[índices locales]...}, de modo que si quiero una lista con los
            índices globales correspondientes a los tres desplazamientos del nodo 200020 y los tres giros del
            200034 debo llamar como:
                nodeDof2idxEng(es.stru, {200020:[1,2,3], 200034:[4,5,6]}).
            
            Para graficar GLG conviene usar el atributo 'u_mdr', auqne podría usarse 'u_raw'. El primero tiene los
            giros expresados como componentes de un vector rotación, mientras que el segundo no tiene los giros,
            pero sí los ángulos de Euler correspondientes a cada instante (lo que sale de CURVAS).
            Por lo tanto, si quiero graficar el desplazamiento vertical y el giro por torsión de la puntera, debo
            usar primero
                nodeDof2idxEng(es.stru, {200035:[3,4]}) -->> [213, 214]
            y luego
                fig_general(es.stru, [[213],[214]], **{'data_type':[['u_mdr'],['u_mdr']]}).
            
            Finalmente, una forma abreviada para hacer listas que tienen elementos repretidos es la siguiente
                3*['q'] -->> ['q', 'q', 'q']
                3*[['q']] -->> [['q'], ['q'], ['q']]
                [2*['q']] -->> [['q', 'q']]
            con lo que la llamada:
                fig_general(es.stru, [[211],[212],[213],[214],[215],[216]], **{'data_type':[['u_mdr'],['u_mdr'],['u_mdr'],['u_mdr'],['u_mdr'],['u_mdr']]})
            puede reducirse a:
                fig_general(es.stru, [[211],[212],[213],[214],[215],[216]], **{'data_type':6*[['u_mdr']]}).
                
            CUIDADO: no sé si fig_general puede graficar FFTs, pero fig_FFT sí, y podemos usar esa
            
        fig_FFT -----
            Necesita argumentos parecido a fig_general, con algunos cambios, a saber:
            1.- Hasta donde entiendo, no puede graficar más de una línea por sistema de ejes
                (esto sería una limitación de la función que calcula y grafica efectivamente la FFT, que es plt-FFT);
                por lo tanto sólo se pueden pasar listas de lista con un elemento, es decir, vale hacer:
                    fig_FFT(es.stru, [[1],[2],[3],[4]], **{'data_type':4*[['q']], 'x_units':'rad/s'}),
                pero no vale hacer:
                    fig_FFT(es.stru, [[1,2,3,4]], **{'data_type':[4*['q']], 'x_units':'rad/s'}).
            2.- Lo probé para los modos, pero no lo probé para los GLG.
                En caso de querer graficar GLG, data_type debe ser UDS (como dice la ayuda de la función).
                Sería bueno que lo pruebes para los GLG y veas si funciona pasarle índices globales o se le puede pasar
                también en forma de diccionario nodo:índices_locales. Si funciona, genial, completás esta ayuda como
                vengo haciendo yo. Si no funciona, me avisás y trata de arreglarlo así queda andando.
            3.- Creo que te conviene agregar 'x_units':'rad/s' en los kwargs, para que sea más fácil comparar,
                pero vos fijate.
    '''
    
    es = rdBin('lin_noAero', **{'subDir_BIN': bin_path})
    
    #   EJEMPLOS DE LAS EXPLICACIONES DE ARRIBA
    # fig_general(es.stru, [[2]], **{'data_type':[['q']]})
    # fig_general(es.stru, [[2],[3]], **{'data_type':[['q'],['q']]})
    # fig_general(es.stru, [[2,3]], **{'data_type':[['q','q']]});
    # fig_general(es.stru, [[2,3],[10],[12]], **{'data_type':[['q','q'],['q'],['q']]});
    # print( nodeDof2idxEng(es.stru, {200020:[1,2,3], 200034:[4,5,6]}) )
    # print( nodeDof2idxEng(es.stru, {200035:[3,4]}) )
    # fig_general(es.stru, [[213],[214]], **{'data_type':[['u_mdr'],['u_mdr']]})
    # fig_general(es.stru, [[211],[212],[213],[214],[215],[216]], **{'data_type':[['u_mdr'],['u_mdr'],['u_mdr'],['u_mdr'],['u_mdr'],['u_mdr']]})
    # fig_general(es.stru, [[211],[212],[213],[214],[215],[216]], **{'data_type':6*[['u_mdr']]})
    # fig_FFT(es.stru, [[1]], **{'data_type':[['q']], 'x_units':'rad/s'})
    # fig_FFT(es.stru, [[1],[2],[3],[4]], **{'data_type':4*[['q']], 'x_units':'rad/s'})
    
    #   CÁLCULO
    # fig_general(es.stru, [[1,2],[3,4],[5,6]], **{'data_type':3*[['q','q']]})
    # fig_general(es.stru, [[1],[2],[3],[4],[5],[6]], **{'data_type':6*[['q']]})
    


