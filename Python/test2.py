import geneontology as go_lib

print("Chargement...")
go = go_lib.load_OBO("Python/data/projet/go-basic.obo")

root_cc = "GO:0005575" # Cellular Component

print(f"\n--- Diagnostic pour {root_cc} ---")
if root_cc in go.nodes:
    print(f"1. Le noeud existe : OUI")
    preds = go.predecessors(root_cc)
    print(f"2. Nombre d'enfants directs (predecessors) : {len(preds)}")
    if len(preds) > 0:
        print(f"   Exemple d'enfant : {preds[0]} -> {go.nodes[preds[0]].get('name')}")
    else:
        print("   ERREUR : Aucun enfant trouvé (C'est pour ça que la profondeur est 0)")
else:
    print("1. Le noeud existe : NON (C'est le problème)")

print("\n--- Recalcul des profondeurs ---")
depths = go_lib.max_depth(go)
print(depths)
go = go_lib.load_OBO("Python/data/projet/go-basic.obo") # Votre chemin

# On regarde l'enfant de la racine CC qui a été trouvé tout à l'heure
child_id = "GO:0032991"

if child_id in go.nodes:
    print(f"Info sur {child_id} :")
    print(go.nodes[child_id])  # On affiche tout le dictionnaire d'attributs
else:
    print(f"{child_id} n'existe pas ??")
