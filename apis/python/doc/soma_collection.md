<a id="tiledbsc.soma_collection"></a>

# tiledbsc.soma\_collection

<a id="tiledbsc.soma_collection.SOMACollection"></a>

## SOMACollection Objects

```python
class SOMACollection(TileDBGroup)
```

Implements a collection of `SOMA` objects.

<a id="tiledbsc.soma_collection.SOMACollection.__init__"></a>

#### \_\_init\_\_

```python
def __init__(uri: str,
             *,
             name: str = "soco",
             soma_options: Optional[SOMAOptions] = None,
             config: Optional[tiledb.Config] = None,
             ctx: Optional[tiledb.Ctx] = None,
             parent: Optional[TileDBGroup] = None)
```

Create a new `SOMACollection` object. The existing group is opened at the

specified `uri` if one is present, otherwise a new group will be created upon ingest.

**Arguments**:

- `uri`: URI of the TileDB group

<a id="tiledbsc.soma_collection.SOMACollection.__repr__"></a>

#### \_\_repr\_\_

```python
def __repr__() -> str
```

Default display of SOMACollection.

<a id="tiledbsc.soma_collection.SOMACollection.__len__"></a>

#### \_\_len\_\_

```python
def __len__() -> int
```

Implements `len(soco)`. Returns the number of elements in the collection.

<a id="tiledbsc.soma_collection.SOMACollection.add"></a>

#### add

```python
def add(soma: SOMA, relative: Optional[bool] = None) -> None
```

Adds a `SOMA` to the `SOMACollection`.

<a id="tiledbsc.soma_collection.SOMACollection.remove"></a>

#### remove

```python
def remove(soma: Union[SOMA, str]) -> None
```

Removes a `SOMA` from the `SOMACollection`, when invoked as `soco.remove("namegoeshere")`.

<a id="tiledbsc.soma_collection.SOMACollection.__delattr__"></a>

#### \_\_delattr\_\_

```python
def __delattr__(matrix_name: str) -> None
```

Removes a `SOMA` from the `SOMACollection`, when invoked as `del soco.namegoeshere`.

<a id="tiledbsc.soma_collection.SOMACollection.__delitem__"></a>

#### \_\_delitem\_\_

```python
def __delitem__(matrix_name: str) -> None
```

Removes a `SOMA` from the `SOMACollection`, when invoked as `del soco["namegoeshere"]`.

<a id="tiledbsc.soma_collection.SOMACollection.keys"></a>

#### keys

```python
def keys() -> Sequence[str]
```

Returns the names of the SOMAs in the collection.

<a id="tiledbsc.soma_collection.SOMACollection.__iter__"></a>

#### \_\_iter\_\_

```python
def __iter__() -> Iterator[SOMA]
```

Implements `for soma in soco: ...`

<a id="tiledbsc.soma_collection.SOMACollection.__contains__"></a>

#### \_\_contains\_\_

```python
def __contains__(name: str) -> bool
```

Implements `name in soco`

<a id="tiledbsc.soma_collection.SOMACollection.__getitem__"></a>

#### \_\_getitem\_\_

```python
def __getitem__(name: str) -> Optional[SOMA]
```

Returns a `SOMA` element at the given name within the group, or `None` if no such
member exists.  Overloads the `[...]` operator.

<a id="tiledbsc.soma_collection.SOMACollection.query"></a>

#### query

```python
def query(*,
          obs_attrs: Optional[Sequence[str]] = None,
          obs_query_string: Optional[str] = None,
          obs_ids: Optional[Ids] = None,
          var_attrs: Optional[Sequence[str]] = None,
          var_query_string: Optional[str] = None,
          var_ids: Optional[Ids] = None) -> Optional[SOMASlice]
```

Subselects the obs, var, and X/data using the specified queries on obs and var,
concatenating across SOMAs in the collection.  Queries use the TileDB-Py `QueryCondition`
API.

If `obs_query_string` is `None`, the `obs` dimension is not filtered and all of `obs` is
used; similiarly for `var`. Return value of `None` indicates an empty slice.  If `obs_ids`
or `var_ids` are not `None`, they are effectively ANDed into the query.  For example, you
can pass in a known list of `obs_ids`, then use `obs_query_string` to further restrict the
query.

If `obs_attrs` or `var_attrs` are unspecified, slices will take all `obs`/`var` attributes
from their source SOMAs; if they are specified, slices will take the specified `obs`/`var`
attributes.  If all SOMAs in the collection have the same `obs`/`var` attributes, then you
needn't specify these; if they don't, you must.

<a id="tiledbsc.soma_collection.SOMACollection.find_unique_obs_values"></a>

#### find\_unique\_obs\_values

```python
def find_unique_obs_values(obs_label: str) -> Set[str]
```

Given an `obs` label such as `cell_type` or `tissue`, returns a set of unique
values for that label among all SOMAs in the collection.

<a id="tiledbsc.soma_collection.SOMACollection.find_unique_var_values"></a>

#### find\_unique\_var\_values

```python
def find_unique_var_values(var_label: str) -> Set[str]
```

Given an `var` label such as `feature_name`, returns a set of unique values for
that label among all SOMAs in the collection.

