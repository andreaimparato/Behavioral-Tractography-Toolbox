function getUnicNameByIndex(index) {
  if (index >= 0 && index < sfereJson.length) {
      return findConnections(sfereJson[index].unic_name);
  } else {
      return "Indice non valido";
  }
}

function findConnections(searchString) {
  const connected = [];

  connection.forEach(item => {
      if (item.src.includes(searchString) || item.trg.includes(searchString)) {
          // Controlla e aggiungi l'indice del nodo src
          if (!connected.includes(item.src)) {
              const srcIndex = sfereJson.findIndex(sfera => sfera.unic_name === item.src);
              if (srcIndex !== -1) {
                  connected.push(srcIndex);
              }
          }
          // Controlla e aggiungi l'indice del nodo trg
          if (!connected.includes(item.trg)) {
              const trgIndex = sfereJson.findIndex(sfera => sfera.unic_name === item.trg);
              if (trgIndex !== -1) {
                  connected.push(trgIndex);
              }
          }
      }
  });

  // Itera sull'array connected e chiama showHideSphere per ogni elemento
  connected.forEach(index => {
      showHideSphere(index);
  });
}






function raccogliNodi() {
    // Inizializza l'array vuoto
    let listaNodi = [];
    
    // Seleziona tutti gli elementi con la classe '.details-item.node-details-item' all'interno di '.drawer.drawer_nodi'
    const nodi = document.querySelectorAll('.drawer.drawer_nodi .details-item.node-details-item');
    
    // Itera sugli elementi selezionati
    nodi.forEach((nodo) => {
      // Cerca il button all'interno dell'elemento .manage che si trova dentro l'elemento nodo
      const button = nodo.querySelector('.manage button');
      
      // Verifica se il button esiste e ha l'attributo 'xx-data-id'
      if (button && button.getAttribute('xx-data-id')) {
        // Ottiene il valore dell'attributo 'xx-data-id' e lo trasforma in numero
        const dataId = Number(button.getAttribute('xx-data-id'));
        
        // Se il valore è un numero valido, lo aggiunge all'array
        if (!isNaN(dataId)) {
          listaNodi.push(dataId);
        }
      }
    });
  
    return listaNodi;
  }
  
  function applicaShowHideSphere() {
    const listaNodi = raccogliNodi();
    // Itera su ogni elemento dell'array listaNodi
    listaNodi.forEach((nodoId) => {
      // Chiama la funzione showHideSphere con il nodoId corrente
      showHideSphere(nodoId);
    });
  }

  function resetShowHideSphere() {
    const listaNodi = raccogliNodi();
    // Itera su ogni elemento dell'array listaNodi
    listaNodi.forEach((nodoId) => {
      // Chiama la funzione showHideSphere con il nodoId corrente
      if(!sfereJson[nodoId].hide)showHideSphere(nodoId);
    });
  }



  // Funzione per iniettare CSS dinamico
function injectCSS() {
  const css = `
    #popupNodi {
      display: none;
      position: fixed;
      top: 50%;
      /* left: 50%; */
      transform: translate(40px, -50%);
      width: 300px;
      background-color: white;
      border: 1px solid #ccc;
      box-shadow: 0 0 10px rgba(0, 0, 0, 0.1);
      z-index: 1000;
      padding: 20px;
      max-height: 70%;
      overflow: scroll;
    }
    #popupNodi ul {
      list-style-type: none;
      padding: 0;
    }
    ul#listaNodi {
      display: flex;
      flex-direction: column;
    }
    #popupNodi li {
      padding: 4px 8px;
      margin-bottom: 0px;
      cursor: pointer;
      border: 1px solid #dddddd;
      border-bottom: 0px;
      order: 1;
    }
    #popupNodi li[data-hidden="false"] {
      color: #0000ff;
      order: 0;
    }
    #popupNodi li:hover {
      color: blue;
    }
    #popupNodi button {
      position: absolute;
      bottom: 10px;
      right: 10px;
      cursor: pointer;
    }
    #overlay {
      display: none;
      position: fixed;
      top: 0;
      left: 0;
      width: 100%;
      height: 100%;
      background: rgba(0, 0, 0, 0.5);
      z-index: 999;
    }
    button#closePopup {
      position: sticky;
      right: 10px;
      left: initial;
      bottom: initial;
      top: 10px;
    }
    div#bottomButton {
      position: absolute;
      bottom: 10px;
      left: 10px;
      z-index: 9999999;
    }
  `;

  const style = document.createElement('style');
  style.textContent = css;
  document.head.appendChild(style);
}

// Funzione per iniettare HTML dinamico per il popup
function injectHTML() {
  const overlay = document.createElement('div');
  overlay.id = 'overlay';
  
  const popup = document.createElement('div');
  popup.id = 'popupNodi';
  
  const closeButton = document.createElement('button');
  closeButton.id = 'closePopup';
  closeButton.textContent = 'X';
  
  const title = document.createElement('h3');
  title.textContent = 'Seleziona un Nodo';
  
  const nodeList = document.createElement('ul');
  nodeList.id = 'listaNodi';
  
  // Aggiunge gli elementi nel popup
  popup.appendChild(closeButton);
  popup.appendChild(title);
  popup.appendChild(nodeList);
  
  // Aggiunge il popup e l'overlay al body
  document.body.appendChild(overlay);
  document.body.appendChild(popup);
  
  // Listener per il pulsante di chiusura del popup
  closeButton.addEventListener('click', chiudiPopupNodi);
}

// Funzione per mostrare il popup con la lista dei nodi
function mostraPopupNodi() {
  const popup = document.getElementById('popupNodi');
  const overlay = document.getElementById('overlay');
  const listaNodiElement = document.getElementById('listaNodi');
  
  // Pulisce eventuali contenuti precedenti
  listaNodiElement.innerHTML = '';

  // Usa la funzione raccogliNodi per ottenere la lista dei nodi
  const listaNodi = raccogliNodi();

  // Aggiunge ciascun nodo alla lista
  listaNodi.forEach((nodoId) => {
    const li = document.createElement('li');
    let layer = sfereJson[nodoId].layer.split('Layer')[1];
    li.textContent = ` ID ${nodoId} ~ ${sfereJson[nodoId].name.replaceAll('_', ' ')} ~ Layer: ${layer}`;
    li.setAttribute('data-id', nodoId);
    li.setAttribute('data-hidden', sfereJson[nodoId].hide);
    li.setAttribute('data-layer', layer);

    // Aggiungi un listener per il clic
    li.addEventListener('click', function () {
      const nodoId = this.getAttribute('data-id');
      getUnicNameByIndex(nodoId);  // Chiama la funzione con xx-data-id
      chiudiPopupNodi();
    });

    listaNodiElement.appendChild(li);
    
  });

  // Mostra il popup e l'overlay
  popup.style.display = 'block';
  overlay.style.display = 'block';
}

// Funzione per chiudere il popup
function chiudiPopupNodi() {
  const popup = document.getElementById('popupNodi');
  const overlay = document.getElementById('overlay');
  
  popup.style.display = 'none';
  overlay.style.display = 'none';
}

// Funzione per inizializzare tutto
function initPopupNodi() {
  injectCSS();
  injectHTML();

  // Listener per aprire il popup, potrebbe essere collegato a un pulsante nella tua interfaccia
  document.getElementById('openPopupButton').addEventListener('click', mostraPopupNodi);
}

// Inizializza il popup quando il documento è pronto
document.addEventListener('DOMContentLoaded', function() {
  initPopupNodi();
});




// Funzione per iniettare il pulsante 'Apri lista nodi' nel DOM
function injectOpenPopupButton() {
  // Crea il pulsante
  const button = document.createElement('button');
  button.id = 'openPopupButton';
  button.textContent = 'Apri lista nodi';

  // Applica gli stili CSS per il posizionamento
  button.style.position = 'absolute';
  button.style.zIndex = '10000';
  button.style.bottom = '20px';
  button.style.left = '20px';

  // Aggiunge il pulsante al body del documento
  //document.body.appendChild(button);

  // Aggiungi il listener per mostrare il popup dei nodi
  //button.addEventListener('click', mostraPopupNodi);

  let htmlTemplate = `
    <div id="bottomButton" class="">
      <button onclick="mostraPopupNodi()" class="lista">apri lista nodi</button>
      <button onclick="resetShowHideSphere()" class="nascondi">nascondi</button>
      <button onclick="applicaShowHideSphere()" class="mostra">mostra</button>
    </div>`;

    $('body').append(htmlTemplate)
}

// Esegui la funzione per iniettare il pulsante quando il documento è pronto
document.addEventListener('DOMContentLoaded', function() {
  injectOpenPopupButton();
});






function aggiornaVisibilitaSfere(sfereJson) {
  // Seleziona tutti gli elementi li all'interno di ul#listaNodi
  const elementi = document.querySelectorAll('#listaNodi li');
  
  // Loop attraverso ogni elemento li
  elementi.forEach(function(elemento) {
      // Recupera il valore dell'attributo data-id
      const dataId = elemento.getAttribute('data-id');
      
      // Controlla se l'oggetto sfereJson contiene un valore per questo data-id
      if (sfereJson[dataId]) {
          // Ottieni il valore di hidden dall'oggetto sfereJson
          const hiddenValue = sfereJson[dataId].hide;
          
          // Assegna il valore di hidden a data-visible dell'elemento li
          elemento.setAttribute('data-hidden', hiddenValue?hiddenValue:'false');
      }
  });
}



// Richiama la funzione
setInterval(function(){aggiornaVisibilitaSfere(sfereJson)}, 1000);



  